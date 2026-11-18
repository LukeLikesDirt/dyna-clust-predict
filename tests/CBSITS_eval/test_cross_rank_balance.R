#!/usr/bin/env Rscript
# test_cross_rank_balance.R — Test the user's cross-rank max_proportion idea:
# when predicting species-level cutoffs for a coarse parent (kingdom, phylum,
# class, ...), currently the only proportion cap applied (R/subset.R's
# nested_prediction_filter, Filter 2) restricts the single dominant SPECIES to
# max_proportion (0.5). Nothing stops one intervening clade (a phylum within
# kingdom Fungi, a class within phylum Ascomycota, etc.) from dominating the
# pool numerically, so the final random downsample to max_sequences can be
# overwhelmingly composed of whichever lineage happens to have the most
# species-identified records in EUKARYOME -- silently turning a "kingdom
# Fungi" cutoff into something close to "Ascomycota's cutoff" instead.
#
# This script compares, for a set of coarse (parent_rank, dataset) pairs that
# are currently pinned at species-level cutoff=1.0:
#   (a) CURRENT sampling: cap dominant species to 50%, then plain random
#       downsample to N -- reproducing R/subset.R::nested_prediction_filter.
#   (b) PROPOSED sampling: recursively cap every intervening rank between
#       parent_rank and species (e.g. phylum -> class -> order -> family ->
#       genus for a kingdom parent) to max 50% *within each branch*, then cap
#       dominant species, then downsample to the same N.
# ...and reports cutoff/F-measure and composition (dominant-clade %,
# singleton fraction) for both, at iddef=2 (the pipeline's current default --
# the iddef question was already tested and resolved separately; this test
# holds it fixed to isolate the max_proportion question).
#
# N is capped at 5000 (not production's 20000) for local tractability -- see
# test_iddef_pinned_species.R for the same tradeoff and its timing rationale.
#
# This is a TEST ONLY. It does not modify R/subset.R or any pipeline script.
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({
  library(data.table)
})

source("R/utils.R")

fasta_all <- "data/full_ITS/eukaryome_ITS_nocomplex.fasta"
class_all <- "data/full_ITS/eukaryome_ITS_nocomplex.classification"
tmp_dir   <- "tmp/cross_rank_balance_test"
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

start_t <- 0.90; end_t <- 1.0; step_t <- 0.001
N_CAP <- 5000L
MAX_PROPORTION <- 0.5
set.seed(42)

# (parent_rank, dataset) pairs to test: every full_ITS species-level cutoff=1
# row at kingdom/phylum/class level (9 groups -- the coarsest, most affected
# by uncontrolled clade dominance), plus two order-level contrasts: one
# fungal (Hypocreales, genuinely species-diverse -- not singleton-pinned even
# under iddef=1) and one plant (Asparagales, singleton-pinned under both
# iddefs -- included to check whether balancing helps or is powerless against
# pure singleton structure).
targets <- data.table(
  parent_rank = c("kingdom", "phylum", "phylum", "class", "class", "class",
                  "class", "class", "class", "order", "order"),
  dataset = c("Fungi", "Ascomycota", "Bryophyta", "Sordariomycetes",
             "Dothideomycetes", "Bryopsida", "Cupressopsida", "Cycadopsida",
             "Gastropoda", "Hypocreales", "Asparagales")
)
# Drop any not actually at higher_rank==parent_rank in the pinned list (the
# 6 "class" rows from the pinned-72 list are the real set; the extra ones
# above duplicate parent_rank so trim to what's actually class-level pinned)
cutoffs <- fread("data/full_ITS/eukaryome.cutoffs.json.txt", sep = "\t", header = TRUE)
pinned_class <- cutoffs[rank == "species" & higher_rank == "class" & `cut-off` == 1, dataset]
targets <- targets[!(parent_rank == "class" & !dataset %in% pinned_class) |
                   parent_rank != "class"]
cat("Testing", nrow(targets), "(parent_rank, dataset) pairs:\n")
print(as.data.frame(targets))

cls <- fread(class_all, sep = "\t", header = TRUE, colClasses = "character")
cls <- cls[is_identified(species)]

# ── Helpers: port of R/subset.R's cap_dominant_group, plus a recursive
#    extension that applies it within every branch at every intervening rank
#    ("class in EACH phylum") ──────────────────────────────────────────────

cap_dominant_group <- function(df, group_col, max_proportion) {
  if (nrow(df) == 0 || max_proportion >= 1) return(df)
  gc <- df[, .N, by = group_col][order(-N)]
  dominant <- gc[[1, group_col]]
  dominant_prop <- gc[[1, "N"]] / nrow(df)
  if (dominant_prop <= max_proportion) return(df)
  non_dom <- df[get(group_col) != dominant]
  dom     <- df[get(group_col) == dominant]
  max_dom <- floor(nrow(non_dom) * max_proportion / (1 - max_proportion))
  dom     <- dom[sample(.N, max_dom)]
  rbind(non_dom, dom)
}

# Bottom-up: recurse into finer ranks FIRST (balancing e.g. genus within each
# order-within-class-within-phylum branch), and only cap the current rank's
# dominance once its branches' sizes are already finalized by that recursion.
# Capping top-down instead (rank r first, then recurse into r's branches)
# looks right but is NOT idempotent: per-branch shrinkage from the finer-rank
# caps applied afterward is uneven across branches (a branch whose own
# children are already balanced shrinks less than one that needed heavy
# trimming), so the coarse-rank proportion set by an early cap silently
# drifts once its branches are later pruned unevenly -- confirmed by a first
# run of this script where "balanced" phylum share came out HIGHER (85%)
# than the unbalanced baseline (60%) for kingdom Fungi, the opposite of the
# intended effect. Capping after recursion (this version) guarantees the cap
# is checked against final branch sizes.
recursive_balance <- function(df, ranks, max_proportion) {
  if (length(ranks) == 0 || nrow(df) == 0) return(df)
  r <- ranks[1]
  df <- df[is_identified(get(r))]
  if (nrow(df) == 0) return(df)
  sub <- split(df, df[[r]])
  balanced <- rbindlist(lapply(sub, recursive_balance, ranks = ranks[-1], max_proportion = max_proportion))
  cap_dominant_group(balanced, r, max_proportion)
}

# ── F-measure sweep (same union-find pattern as test_iddef_pinned_species.R) ─

sweep_f <- function(sim, class_of, class_size, n) {
  parent <- seq_len(n); rnk <- integer(n); csize <- rep(1L, n)
  n_classes <- length(class_size)
  root_cc <- lapply(seq_len(n), function(i) {
    v <- integer(n_classes); if (class_of[i] > 0L) v[class_of[i]] <- 1L; v
  })
  is_root <- rep(TRUE, n)
  find <- function(x) { while (parent[x] != x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }; x }
  union <- function(a, b) {
    ra <- find(a); rb <- find(b); if (ra == rb) return(invisible())
    if (rnk[ra] < rnk[rb]) { t <- ra; ra <- rb; rb <- t }
    parent[rb] <<- ra
    if (rnk[ra] == rnk[rb]) rnk[ra] <<- rnk[ra] + 1L
    csize[ra] <<- csize[ra] + csize[rb]
    root_cc[[ra]] <<- root_cc[[ra]] + root_cc[[rb]]
    is_root[rb] <<- FALSE
  }
  fm <- function() {
    active <- which(is_root); tot <- 0
    for (ci in seq_len(n_classes)) {
      if (class_size[ci] == 0L) next
      best <- 0
      for (r in active) {
        cnt <- root_cc[[r]][ci]
        if (cnt > 0L) { d <- 2 * cnt / (class_size[ci] + csize[r]); if (d > best) best <- d }
      }
      tot <- tot + class_size[ci] * best
    }
    tot / sum(class_size)
  }
  if (nrow(sim) == 0) return(list(cut = NA, f = NA))
  ord <- order(-sim$score); ii <- sim$i[ord]; jj <- sim$j[ord]; ss <- sim$score[ord]
  ptr <- 1L; ne <- length(ss)
  ts <- rev(seq(start_t, end_t, by = step_t)); fs <- numeric(length(ts))
  for (k in seq_along(ts)) {
    t <- ts[k]
    while (ptr <= ne && ss[ptr] >= t - 1e-9) { union(ii[ptr], jj[ptr]); ptr <- ptr + 1L }
    fs[k] <- fm()
  }
  ts <- rev(ts); fs <- rev(fs)
  best <- max(fs); tied <- which(fs == best)
  list(cut = ts[tied[ceiling(length(tied) / 2)]], f = round(best, 4))
}

run_vsearch_f <- function(ids, sp_map, safe_name, n_threads) {
  fa_file  <- file.path(tmp_dir, paste0(safe_name, ".fasta"))
  id_file  <- file.path(tmp_dir, paste0(safe_name, ".ids"))
  sim_file <- file.path(tmp_dir, paste0(safe_name, ".sim"))
  writeLines(ids, id_file)
  system(sprintf("seqkit grep -f '%s' '%s' > '%s' 2>/dev/null", id_file, fasta_all, fa_file))
  system(sprintf(paste("vsearch --allpairs_global '%s' --acceptall --userout '%s'",
                       "--userfields query+target+id --iddef 2 --threads %d"),
                 fa_file, sim_file, n_threads),
         ignore.stdout = TRUE, ignore.stderr = TRUE)
  res <- list(cut = NA, f = NA)
  if (file.exists(sim_file) && file.size(sim_file) > 0) {
    raw <- fread(sim_file, header = FALSE, col.names = c("q", "t", "pid"),
                colClasses = c("character", "character", "numeric"))
    raw[, score := pid / 100]
    idx <- setNames(seq_along(ids), ids)
    sp <- sp_map[ids]
    classes <- unique(sp[!is.na(sp)])
    class_of <- match(sp, classes); class_of[is.na(class_of)] <- 0L
    class_size <- tabulate(class_of[class_of > 0], nbins = length(classes))
    res <- sweep_f(data.table(i = idx[raw$q], j = idx[raw$t], score = raw$score),
                   class_of, class_size, length(ids))
  }
  unlink(c(fa_file, id_file, sim_file))
  res
}

# ── Composition summary: dominant share at each rank, singleton fraction ────

composition_summary <- function(df, ranks) {
  parts <- sapply(ranks, function(r) {
    df2 <- df[is_identified(get(r))]
    if (nrow(df2) == 0) return(NA_real_)
    max(df2[, .N, by = r]$N) / nrow(df2)
  })
  n_spp <- length(unique(df$species))
  sprintf("%s | n_spp/n_seq=%d/%d (%.1f%% singleton-rate proxy)",
         paste(sprintf("%s=%.0f%%", ranks, parts * 100), collapse = ", "),
         n_spp, nrow(df), 100 * (1 - (nrow(df) - n_spp) / nrow(df)))
}

# ── Process each target ──────────────────────────────────────────────────────

out <- list()
for (i in seq_len(nrow(targets))) {
  pr <- targets$parent_rank[i]; ds <- targets$dataset[i]
  intervening <- setdiff(parent_ranks_map[["species"]], pr)
  # parent_ranks_map$species is finest->coarsest genus..kingdom; keep only
  # ranks strictly between pr and species (coarser than genus cut where pr
  # itself sits), ordered coarse->fine to match "phylum, then class within
  # phylum, then order within class..." recursion order.
  full_order <- c("phylum", "class", "order", "family", "genus")
  start_idx <- if (pr == "kingdom") 1L else which(full_order == pr) + 1L
  ranks_between <- if (pr %in% full_order) full_order[start_idx:length(full_order)] else full_order
  if (pr == "kingdom") ranks_between <- full_order

  pool <- cls[get(pr) == ds]
  if (nrow(pool) < 30) { cat("skip (too few):", pr, ds, "\n"); next }
  sp_map <- setNames(pool$species, pool$id)

  cat(sprintf("\n--- %s / %s (pool n=%d, balancing ranks: %s) ---\n",
              pr, ds, nrow(pool), paste(ranks_between, collapse = " > ")))

  # (a) CURRENT: cap dominant species only, then random downsample
  cur <- cap_dominant_group(pool, "species", MAX_PROPORTION)
  if (nrow(cur) > N_CAP) cur <- cur[sample(.N, N_CAP)]

  # (b) PROPOSED: recursively balance intervening ranks, then cap species,
  #     then downsample to the same N
  bal <- recursive_balance(pool, ranks_between, MAX_PROPORTION)
  bal <- cap_dominant_group(bal, "species", MAX_PROPORTION)
  if (nrow(bal) > N_CAP) bal <- bal[sample(.N, N_CAP)]

  cat("  current : ", composition_summary(cur, ranks_between), "\n")
  cat("  balanced: ", composition_summary(bal, ranks_between), "\n")

  n_threads <- if (max(nrow(cur), nrow(bal)) > 500) 10L else 4L
  safe <- gsub("[^A-Za-z0-9_]", "_", paste(pr, ds, sep = "_"))

  t0 <- Sys.time()
  r_cur <- run_vsearch_f(cur$id, sp_map, paste0(safe, "_current"), n_threads)
  r_bal <- run_vsearch_f(bal$id, sp_map, paste0(safe, "_balanced"), n_threads)
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)

  cat(sprintf("  [%d/%d] (%.1fs) current: cut=%s F=%s | balanced: cut=%s F=%s\n",
              i, nrow(targets), elapsed,
              ifelse(is.na(r_cur$cut), "NA", sprintf("%.3f", r_cur$cut)),
              ifelse(is.na(r_cur$f), "NA", sprintf("%.4f", r_cur$f)),
              ifelse(is.na(r_bal$cut), "NA", sprintf("%.3f", r_bal$cut)),
              ifelse(is.na(r_bal$f), "NA", sprintf("%.4f", r_bal$f))))

  out[[length(out) + 1]] <- data.table(
    parent_rank = pr, dataset = ds,
    n_current = nrow(cur), n_spp_current = length(unique(cur$species)),
    cut_current = r_cur$cut, f_current = r_cur$f,
    n_balanced = nrow(bal), n_spp_balanced = length(unique(bal$species)),
    cut_balanced = r_bal$cut, f_balanced = r_bal$f,
    elapsed_sec = elapsed
  )
  fwrite(rbindlist(out), "tests/CBSITS_eval/results/cross_rank_balance.tsv", sep = "\t")
}

res <- rbindlist(out)
cat("\n=== SUMMARY ===\n")
print(as.data.frame(res), row.names = FALSE)
cat("\nWrote tests/CBSITS_eval/results/cross_rank_balance.tsv\n")
