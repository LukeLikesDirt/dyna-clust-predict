#!/usr/bin/env Rscript
# test_singleton_70_30.R — Gentler version of the max-singleton-fraction
# constraint tested in test_singleton_constraints.R: instead of capping
# singleton species at 50% of the species count, cap at 70% (a 70/30
# singleton:multiseq ratio). Motivation: at 50%, kingdom Fungi's optimum
# landed at cutoff~0.95, but real fungal ITS is expected closer to ~0.98
# (see tests/CBSITS_eval/fungi_curve.tsv / the published curve artifact --
# the F-measure curve is on a flat plateau from ~0.94-0.965 and only starts
# declining meaningfully above ~0.97, so a gentler singleton constraint that
# keeps more real data might land closer to that plateau's upper shoulder).
#
# Same 14 (parent_rank, dataset) targets as test_singleton_constraints.R, so
# results merge directly into one comparison table. Reuses (a) current
# baseline from cross_rank_balance.tsv and (b) the already-computed 50%
# singleton-cap results from singleton_constraints.tsv -- only the new 70%
# arm is computed here.
#
# This is a TEST ONLY. It does not modify R/predict.R, R/subset.R, or any
# pipeline script.
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({
  library(data.table)
})

source("R/utils.R")

fasta_all <- "data/full_ITS/eukaryome_ITS_nocomplex.fasta"
class_all <- "data/full_ITS/eukaryome_ITS_nocomplex.classification"
tmp_dir   <- "tmp/singleton_70_30_test"
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

start_t <- 0.90; end_t <- 1.0; step_t <- 0.001
N_CAP <- 5000L
MAX_PROPORTION <- 0.5
MAX_SINGLETON_FRAC <- 0.7
set.seed(42)

targets <- data.table(
  parent_rank = c("kingdom", "phylum", "phylum", "class", "class", "class",
                  "class", "class", "class", "order", "order",
                  "genus", "genus", "genus"),
  dataset = c("Fungi", "Ascomycota", "Bryophyta", "Sordariomycetes",
             "Dothideomycetes", "Bryopsida", "Cupressopsida", "Cycadopsida",
             "Gastropoda", "Hypocreales", "Asparagales",
             "Aspergillus", "Fusarium", "Colletotrichum")
)
cat("Testing", nrow(targets), "(parent_rank, dataset) pairs at 70/30 singleton ratio:\n")
print(as.data.frame(targets))

cls <- fread(class_all, sep = "\t", header = TRUE, colClasses = "character")
cls <- cls[is_identified(species)]

baseline_current <- fread("tests/CBSITS_eval/results/cross_rank_balance.tsv", sep = "\t", header = TRUE)
baseline_current <- baseline_current[, .(parent_rank, dataset, n_current, n_spp_current, cut_current, f_current)]
baseline_50 <- fread("tests/CBSITS_eval/results/singleton_constraints.tsv", sep = "\t", header = TRUE)
baseline_50 <- baseline_50[, .(parent_rank, dataset, n_max50, n_spp_max50, pct_singleton_max50, cut_max50, f_max50)]

# ── Pool-construction helpers (identical to test_singleton_constraints.R) ───

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

cap_singleton_fraction <- function(df, max_singleton_frac) {
  spp_counts <- df[, .N, by = species]
  singleton_spp <- spp_counts[N == 1, species]
  multi_spp     <- spp_counts[N >= 2, species]
  n_multi <- length(multi_spp); n_single <- length(singleton_spp)
  if (n_single == 0 || n_multi == 0) return(df)
  max_single_kept <- floor(max_singleton_frac / (1 - max_singleton_frac) * n_multi)
  if (n_single <= max_single_kept) return(df)
  keep_singleton <- sample(singleton_spp, max_single_kept)
  df[species %in% multi_spp | species %in% keep_singleton]
}

downsample_by_species <- function(df, n_cap) {
  if (nrow(df) <= n_cap) return(df)
  counts <- df[, .N, by = species]
  spp <- sample(counts$species)
  setkey(counts, species)
  cum <- 0L; keep <- character(0)
  for (s in spp) {
    ns <- counts[s, N]
    if (cum + ns > n_cap && length(keep) > 0) next
    keep <- c(keep, s); cum <- cum + ns
    if (cum >= n_cap) break
  }
  df[species %in% keep]
}

# ── F-measure sweep ───────────────────────────────────────────────────────

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

singleton_pct <- function(df) {
  spp_counts <- df[, .N, by = species]
  100 * sum(spp_counts$N == 1) / nrow(spp_counts)
}

# ── Process each target ──────────────────────────────────────────────────────

out <- list()
for (i in seq_len(nrow(targets))) {
  pr <- targets$parent_rank[i]; ds <- targets$dataset[i]
  pool <- cls[get(pr) == ds]
  if (nrow(pool) < 30) { cat("skip (too few):", pr, ds, "\n"); next }
  sp_map <- setNames(pool$species, pool$id)

  cat(sprintf("\n--- %s / %s (raw pool n=%d, %.1f%% singleton species) ---\n",
              pr, ds, nrow(pool), singleton_pct(pool)))

  base_pool <- cap_dominant_group(pool, "species", MAX_PROPORTION)
  d_pool <- cap_singleton_fraction(base_pool, MAX_SINGLETON_FRAC)
  d_pool <- downsample_by_species(d_pool, N_CAP)

  cat(sprintf("  max70pct_singleton: n=%d, %.1f%% singleton spp, %d spp\n",
              nrow(d_pool), singleton_pct(d_pool), length(unique(d_pool$species))))

  n_threads <- if (nrow(d_pool) > 500) 10L else 4L
  safe <- gsub("[^A-Za-z0-9_]", "_", paste(pr, ds, sep = "_"))

  t0 <- Sys.time()
  r_d <- if (nrow(d_pool) >= 30) run_vsearch_f(d_pool$id, sp_map, paste0(safe, "_d"), n_threads) else list(cut = NA, f = NA)
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)

  cat(sprintf("  [%d/%d] (%.1fs) max70pct_singleton: cut=%s F=%s\n",
              i, nrow(targets), elapsed,
              ifelse(is.na(r_d$cut), "NA", sprintf("%.3f", r_d$cut)), ifelse(is.na(r_d$f), "NA", sprintf("%.4f", r_d$f))))

  out[[length(out) + 1]] <- data.table(
    parent_rank = pr, dataset = ds,
    n_max70 = nrow(d_pool), n_spp_max70 = length(unique(d_pool$species)),
    pct_singleton_max70 = round(singleton_pct(d_pool), 1),
    cut_max70 = r_d$cut, f_max70 = r_d$f,
    elapsed_sec = elapsed
  )
  fwrite(rbindlist(out), "tests/CBSITS_eval/results/singleton_70_30.tsv", sep = "\t")
}

res70 <- rbindlist(out)
merged <- baseline_current[res70, on = c("parent_rank", "dataset")]
merged <- baseline_50[merged, on = c("parent_rank", "dataset")]
setcolorder(merged, c("parent_rank", "dataset",
                      "n_current", "n_spp_current", "cut_current", "f_current",
                      "n_max70", "n_spp_max70", "pct_singleton_max70", "cut_max70", "f_max70",
                      "n_max50", "n_spp_max50", "pct_singleton_max50", "cut_max50", "f_max50"))

cat("\n=== SUMMARY: current vs 70/30 vs 50/50 singleton cap ===\n")
print(as.data.frame(merged), row.names = FALSE)
fwrite(merged, "tests/CBSITS_eval/results/singleton_70_30_merged.tsv", sep = "\t")
cat("\nWrote tests/CBSITS_eval/results/singleton_70_30.tsv and singleton_70_30_merged.tsv\n")
