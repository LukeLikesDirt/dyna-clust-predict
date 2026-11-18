#!/usr/bin/env Rscript
# test_fungi_grid.R — Grid test on kingdom Fungi (the case that exposed the
# downsampling-manufactures-fake-singletons finding): cross max sequences
# per species x max singleton-species proportion, both applied on top of
# species-aware (whole-species) downsampling instead of production's
# row-level slice_sample, which was shown to inflate the sample's apparent
# singleton rate from 50.3% (raw) to 89.0% (current 5000-seq sample) by
# fragmenting real multi-sequence species down to 1 surviving sequence each.
#
# Rationale for the max-seqs-per-species axis: a species with many sequences
# and real intraspecific diversity can't all be pairwise-identical at a
# threshold near 1.0, so it has to fragment into sub-clusters there --
# lowering its own Dice contribution at high thresholds and creating
# downward pressure on the optimum (unlike a singleton, which only ever
# rewards staying at 1.0 and never penalizes it). But that pull is only
# useful when it reflects real biology, not database oversampling: 3.3% of
# Fungi species (1,414 of 42,556) have >=10 sequences yet account for 26.1%
# of all sequences (Thanatephorus cucumeris alone has 452), almost certainly
# reflecting medical/agricultural sequencing interest rather than uniform
# biological sampling. Capping sequences-per-species bounds how much any one
# such species -- oversampled for non-biological reasons -- can dominate the
# sequence-count-weighted F sum, independent of the singleton-proportion axis.
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
tmp_dir   <- "tmp/fungi_grid_test"
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

start_t <- 0.90; end_t <- 1.0; step_t <- 0.001
N_CAP <- 5000L
MAX_PROPORTION <- 0.5   # existing species-dominance cap, held fixed
set.seed(42)

MAX_SEQS_GRID      <- c(3, 5, 10, Inf)
MAX_SINGLETON_GRID <- c(0.3, 0.5, 0.7, 1.0)   # 1.0 = no cap

cls <- fread(class_all, sep = "\t", header = TRUE, colClasses = "character")
cls <- cls[is_identified(species)]
pool <- cls[kingdom == "Fungi"]
sp_map <- setNames(pool$species, pool$id)
cat(sprintf("Fungi raw pool: n=%d sequences, %d species\n", nrow(pool), length(unique(pool$species))))

# ── Pool-construction helpers ────────────────────────────────────────────────

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

# Cap how many sequences any single species can contribute: species over the
# cap are randomly thinned down to it (kept, not dropped -- unlike the
# singleton-proportion cap, this never removes a species entirely).
cap_seqs_per_species <- function(df, max_seqs) {
  if (!is.finite(max_seqs)) return(df)
  counts <- df[, .N, by = species]
  over <- counts[N > max_seqs, species]
  if (length(over) == 0) return(df)
  under_df <- df[!species %in% over]
  over_df  <- df[species %in% over]
  thinned <- over_df[, .SD[sample(.N, max_seqs)], by = species]
  rbind(under_df, thinned)
}

cap_singleton_fraction <- function(df, max_singleton_frac) {
  if (max_singleton_frac >= 1) return(df)
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

base_pool <- cap_dominant_group(pool, "species", MAX_PROPORTION)

# ── Grid ──────────────────────────────────────────────────────────────────

out <- list()
combo_i <- 0
n_combos <- length(MAX_SEQS_GRID) * length(MAX_SINGLETON_GRID)
for (max_seqs in MAX_SEQS_GRID) {
  for (max_single in MAX_SINGLETON_GRID) {
    combo_i <- combo_i + 1
    g_pool <- cap_seqs_per_species(base_pool, max_seqs)
    g_pool <- cap_singleton_fraction(g_pool, max_single)
    g_pool <- downsample_by_species(g_pool, N_CAP)

    n_threads <- if (nrow(g_pool) > 500) 10L else 4L
    safe <- sprintf("seq%s_single%s", ifelse(is.finite(max_seqs), max_seqs, "Inf"),
                    gsub("\\.", "", sprintf("%.1f", max_single)))

    t0 <- Sys.time()
    r <- if (nrow(g_pool) >= 30) run_vsearch_f(g_pool$id, sp_map, safe, n_threads) else list(cut = NA, f = NA)
    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)

    pct_s <- singleton_pct(g_pool)
    cat(sprintf("[%2d/%2d] max_seqs=%-4s max_single=%.1f | n=%5d n_spp=%5d singleton=%.1f%% | cut=%s F=%s (%.1fs)\n",
                combo_i, n_combos, ifelse(is.finite(max_seqs), max_seqs, "Inf"), max_single,
                nrow(g_pool), length(unique(g_pool$species)), pct_s,
                ifelse(is.na(r$cut), "NA", sprintf("%.3f", r$cut)),
                ifelse(is.na(r$f), "NA", sprintf("%.4f", r$f)), elapsed))

    out[[length(out) + 1]] <- data.table(
      max_seqs_per_species = ifelse(is.finite(max_seqs), max_seqs, NA_real_),
      max_singleton_prop = max_single,
      n = nrow(g_pool), n_spp = length(unique(g_pool$species)),
      pct_singleton = round(pct_s, 1),
      cutoff = r$cut, f_measure = r$f, elapsed_sec = elapsed
    )
    fwrite(rbindlist(out), "tests/CBSITS_eval/results/fungi_grid.tsv", sep = "\t")
  }
}

res <- rbindlist(out)
cat("\n=== FULL GRID ===\n")
print(as.data.frame(res), row.names = FALSE)

cat("\n=== cutoff, as a matrix (rows=max_seqs_per_species, cols=max_singleton_prop) ===\n")
print(dcast(res, max_seqs_per_species ~ max_singleton_prop, value.var = "cutoff"))
cat("\n=== F-measure, as a matrix ===\n")
print(dcast(res, max_seqs_per_species ~ max_singleton_prop, value.var = "f_measure"))

cat("\nWrote tests/CBSITS_eval/results/fungi_grid.tsv\n")
