#!/usr/bin/env Rscript
# curve_fungi_singleton.R — Re-derive the kingdom Fungi pools exactly as
# test_singleton_constraints.R built them (same seed=42, same operation
# order, Fungi is the first target processed so no prior RNG draws differ),
# but this time keep the FULL per-threshold F-measure curve instead of just
# the argmax, so the peak shape and the cost of moving off it can be
# inspected directly -- requested before deciding on a 70/30 singleton-ratio
# test, since 0.95 (what max50pct_singleton/min2seq landed on) reads too low
# against the ~0.98 expected for real fungal ITS.
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({
  library(data.table)
})

source("R/utils.R")

fasta_all <- "data/full_ITS/eukaryome_ITS_nocomplex.fasta"
class_all <- "data/full_ITS/eukaryome_ITS_nocomplex.classification"
tmp_dir   <- "tmp/curve_fungi_test"
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

start_t <- 0.90; end_t <- 1.0; step_t <- 0.001
N_CAP <- 5000L
MAX_PROPORTION <- 0.5
MAX_SINGLETON_FRAC <- 0.5
set.seed(42)

cls <- fread(class_all, sep = "\t", header = TRUE, colClasses = "character")
cls <- cls[is_identified(species)]

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
drop_singletons <- function(df) {
  spp_counts <- df[, .N, by = species]
  multi_spp <- spp_counts[N >= 2, species]
  df[species %in% multi_spp]
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

pool <- cls[kingdom == "Fungi"]
base_pool <- cap_dominant_group(pool, "species", MAX_PROPORTION)
b_pool <- cap_singleton_fraction(base_pool, MAX_SINGLETON_FRAC)
b_pool <- downsample_by_species(b_pool, N_CAP)
c_pool <- drop_singletons(base_pool)
c_pool <- downsample_by_species(c_pool, N_CAP)

spc_b <- b_pool[, .N, by = species]
spc_c <- c_pool[, .N, by = species]
cat(sprintf("max50pct_singleton repro: n=%d, %.1f%% singleton, %d spp  (expect n=5000, 50.6%%, 1687 spp)\n",
            nrow(b_pool), 100 * mean(spc_b$N == 1), nrow(spc_b)))
cat(sprintf("min2seq repro           : n=%d, %.1f%% singleton, %d spp  (expect n=4999, 0.0%%, 1132 spp)\n",
            nrow(c_pool), 100 * mean(spc_c$N == 1), nrow(spc_c)))

# ── Full-trace F sweep (same union-find pattern, but keeps every threshold) ──

sweep_f_full <- function(sim, class_of, class_size, n) {
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
  if (nrow(sim) == 0) return(data.table(threshold = numeric(0), f = numeric(0)))
  ord <- order(-sim$score); ii <- sim$i[ord]; jj <- sim$j[ord]; ss <- sim$score[ord]
  ptr <- 1L; ne <- length(ss)
  ts <- rev(seq(start_t, end_t, by = step_t)); fs <- numeric(length(ts))
  for (k in seq_along(ts)) {
    t <- ts[k]
    while (ptr <= ne && ss[ptr] >= t - 1e-9) { union(ii[ptr], jj[ptr]); ptr <- ptr + 1L }
    fs[k] <- fm()
  }
  data.table(threshold = rev(ts), f = round(rev(fs), 4))
}

run_curve <- function(pool_df, safe_name, n_threads) {
  ids <- pool_df$id
  sp_map <- setNames(pool_df$species, pool_df$id)
  fa_file  <- file.path(tmp_dir, paste0(safe_name, ".fasta"))
  id_file  <- file.path(tmp_dir, paste0(safe_name, ".ids"))
  sim_file <- file.path(tmp_dir, paste0(safe_name, ".sim"))
  writeLines(ids, id_file)
  system(sprintf("seqkit grep -f '%s' '%s' > '%s' 2>/dev/null", id_file, fasta_all, fa_file))
  system(sprintf(paste("vsearch --allpairs_global '%s' --acceptall --userout '%s'",
                       "--userfields query+target+id --iddef 2 --threads %d"),
                 fa_file, sim_file, n_threads),
         ignore.stdout = TRUE, ignore.stderr = TRUE)
  raw <- fread(sim_file, header = FALSE, col.names = c("q", "t", "pid"),
              colClasses = c("character", "character", "numeric"))
  raw[, score := pid / 100]
  idx <- setNames(seq_along(ids), ids)
  sp <- sp_map[ids]
  classes <- unique(sp[!is.na(sp)])
  class_of <- match(sp, classes); class_of[is.na(class_of)] <- 0L
  class_size <- tabulate(class_of[class_of > 0], nbins = length(classes))
  curve <- sweep_f_full(data.table(i = idx[raw$q], j = idx[raw$t], score = raw$score),
                        class_of, class_size, length(ids))
  unlink(c(fa_file, id_file, sim_file))
  curve
}

cat("\nRunning vsearch for max50pct_singleton pool...\n")
curve_b <- run_curve(b_pool, "fungi_max50", 12L)
cat("Running vsearch for min2seq pool...\n")
curve_c <- run_curve(c_pool, "fungi_min2seq", 12L)

curve_b[, condition := "max50pct_singleton"]
curve_c[, condition := "min2seq"]
curves <- rbind(curve_b, curve_c)
fwrite(curves, "tests/CBSITS_eval/results/fungi_curve.tsv", sep = "\t")

print_curve <- function(curve, label) {
  cat(sprintf("\n=== %s ===\n", label))
  peak_f <- max(curve$f)
  peak_t <- curve[f == peak_f, threshold]
  cat(sprintf("Peak F = %.4f at threshold(s): %s\n", peak_f, paste(sprintf("%.3f", peak_t), collapse = ", ")))
  show <- curve[threshold %in% seq(0.90, 1.00, by = 0.005)]
  setorder(show, -threshold)
  show[, delta_from_peak := round(f - peak_f, 4)]
  print(as.data.frame(show), row.names = FALSE)
}

print_curve(curve_b, "max50pct_singleton: F vs threshold (0.005 steps)")
print_curve(curve_c, "min2seq: F vs threshold (0.005 steps)")

cat("\nWrote tests/CBSITS_eval/results/fungi_curve.tsv (full 0.001-step traces, both conditions)\n")
