#!/usr/bin/env Rscript
# test_random_vs_species_selection.R — kingdom Fungi, --iddef 2, subset.R-style
# limits (min_subgroups 10, min_sequences 30, max_proportion 0.7, max_sequences
# N), NO per-species cap and NO singleton-proportion cap. Compares how the
# sample is drawn when the pool exceeds max_sequences:
#   random_rows   : plain random sequence (row-level) sample from the full pool
#   species_aware : whole species kept/dropped (production's downsample_by_group)
#
# Faster F-measure sweep than test_fungi_grid.R: vsearch only reports pairs
# >= 0.90 identity (--id), and each threshold's single-linkage clusters come
# from igraph components with a vectorised Dice calculation. Semantics match
# sweep_f() in test_fungi_grid.R (validated in validate_fast_sweep.R).
#
# This is a TEST ONLY. It does not modify R/predict.R, R/subset.R, or any
# pipeline script.
#
# Usage (from project root): Rscript tests/CBSITS_eval/test_random_vs_species_selection.R <N> <seeds,comma> 
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({ library(data.table); library(igraph) })
source("R/utils.R")

args <- commandArgs(trailingOnly = TRUE)
N_CAP <- as.integer(if (length(args) >= 1) args[1] else 5000)
SEEDS <- as.integer(strsplit(if (length(args) >= 2) args[2] else "42,43,44", ",")[[1]])
ARMS <- strsplit(if (length(args) >= 3) args[3] else "random_rows,species_aware", ",")[[1]]
ARM_TAG <- if (length(args) >= 3) paste0("_", gsub("[^a-z0-9]", "", args[3])) else ""
MAX_SINGLETON_PROP <- 0.7   # used only by the random_rows_cap70 arm
MIN_SUBGROUPS <- 10L; MIN_SEQUENCES <- 30L; MAX_PROPORTION <- 0.7
IDDEF <- 2L; N_THREADS <- 10L
start_t <- 0.90; end_t <- 1.0; step_t <- 0.001

fasta_all <- "data/full_ITS/eukaryome_ITS_nocomplex.fasta"
class_all <- "data/full_ITS/eukaryome_ITS_nocomplex.classification"
tmp_dir <- file.path("tmp", sprintf("random_vs_species_N%d", N_CAP)); dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- sprintf("tests/CBSITS_eval/results/random_vs_species_selection_N%d%s.tsv", N_CAP, ARM_TAG)

# Drop randomly chosen singleton species (never multi-sequence ones) until they
# are at most max_frac of the species count. Applied AFTER row sampling in the
# random_rows_cap70 arm; the sample therefore ends up smaller than N_CAP.
cap_singleton_fraction <- function(df, max_frac) {
  if (max_frac >= 1) return(df)
  cnt <- df[, .N, by = species]
  single <- cnt[N == 1, species]; multi <- cnt[N >= 2, species]
  if (length(single) == 0 || length(multi) == 0) return(df)
  max_keep <- floor(max_frac / (1 - max_frac) * length(multi))
  if (length(single) <= max_keep) return(df)
  keep <- sample(single, max_keep)
  df[species %in% multi | species %in% keep]
}

cap_dominant_group <- function(df, group_col, max_proportion) {
  if (nrow(df) == 0 || max_proportion >= 1) return(df)
  gc <- df[, .N, by = group_col][order(-N)]
  dominant <- gc[[1, group_col]]
  if (gc[[1, "N"]] / nrow(df) <= max_proportion) return(df)
  non_dom <- df[get(group_col) != dominant]; dom <- df[get(group_col) == dominant]
  rbind(non_dom, dom[sample(.N, floor(nrow(non_dom) * max_proportion / (1 - max_proportion)))])
}
downsample_by_species <- function(df, n_cap) {
  if (nrow(df) <= n_cap) return(df)
  counts <- df[, .N, by = species]; spp <- sample(counts$species); setkey(counts, species)
  cum <- 0L; keep <- character(0)
  for (s in spp) { ns <- counts[s, N]; if (cum + ns > n_cap && length(keep) > 0) next
    keep <- c(keep, s); cum <- cum + ns; if (cum >= n_cap) break }
  df[species %in% keep]
}

# Fast sweep: edges = data.table(i, j, score) with score >= start_t only
fast_sweep <- function(edges, class_of, class_size, n) {
  setorder(edges, -score)
  ts <- seq(start_t, end_t, by = step_t); fs <- numeric(length(ts))
  csz <- class_size; tot <- sum(csz)
  for (k in seq_along(ts)) {
    e <- edges[score >= ts[k] - 1e-9]
    memb <- if (nrow(e) == 0) seq_len(n) else
      components(graph_from_edgelist(cbind(e$i, e$j), directed = FALSE) |> (\(g) add_vertices(g, n - vcount(g)))())$membership
    dt <- data.table(root = memb, cls = class_of)[cls > 0L, .N, by = .(root, cls)]
    rsz <- tabulate(memb, nbins = max(memb))
    dt[, d := 2 * N / (csz[cls] + rsz[root])]
    best <- dt[, .(b = max(d)), by = cls]
    fs[k] <- sum(csz[best$cls] * best$b) / tot
  }
  b <- max(fs); tied <- which(fs == b)
  list(cut = ts[tied[ceiling(length(tied) / 2)]], f = round(b, 4), f_at_1 = round(fs[length(fs)], 4))
}

run_one <- function(ids, sp_map, name) {
  fa <- file.path(tmp_dir, paste0(name, ".fasta")); idf <- file.path(tmp_dir, paste0(name, ".ids"))
  sim <- file.path(tmp_dir, paste0(name, ".sim")); writeLines(ids, idf)
  system(sprintf("seqkit grep -f '%s' '%s' > '%s' 2>/dev/null", idf, fasta_all, fa))
  system(sprintf("vsearch --allpairs_global '%s' --id %.2f --userout '%s' --userfields query+target+id --iddef %d --threads %d",
                 fa, start_t, sim, IDDEF, N_THREADS), ignore.stdout = TRUE, ignore.stderr = TRUE)
  raw <- fread(sim, header = FALSE, col.names = c("q", "t", "pid"), colClasses = c("character", "character", "numeric"))
  idx <- setNames(seq_along(ids), ids); sp <- sp_map[ids]
  classes <- unique(sp); class_of <- match(sp, classes)
  class_size <- tabulate(class_of, nbins = length(classes))
  r <- fast_sweep(data.table(i = idx[raw$q], j = idx[raw$t], score = raw$pid / 100), class_of, class_size, length(ids))
  r$n_edges <- nrow(raw); unlink(c(fa, idf, sim)); r
}

if (sys.nframe() == 0) {
  cls <- fread(class_all, sep = "\t", header = TRUE, colClasses = "character")
  cls <- cls[is_identified(species)]; pool <- cls[kingdom == "Fungi"]
  sp_map <- setNames(pool$species, pool$id)
  raw_n <- pool[, .N, by = species]; multi_in_raw <- raw_n[N >= 2, species]
  cat(sprintf("Fungi pool: n=%d, %d species, %.1f%% singleton by species\n",
              nrow(pool), nrow(raw_n), 100 * mean(raw_n$N == 1)))
  out <- list()
  for (seed in SEEDS) for (arm in ARMS) {
    set.seed(seed); base <- cap_dominant_group(pool, "species", MAX_PROPORTION)
    set.seed(seed + 1000L)
    s <- if (nrow(base) <= N_CAP) base else if (arm %in% c("random_rows", "random_rows_cap70")) base[sample(.N, N_CAP)] else downsample_by_species(base, N_CAP)
    if (arm == "random_rows_cap70") s <- cap_singleton_fraction(s, MAX_SINGLETON_PROP)
    sc <- s[, .N, by = species]
    stopifnot(nrow(sc) >= MIN_SUBGROUPS, nrow(s) >= MIN_SEQUENCES)
    fake <- 100 * sum(sc$N == 1 & sc$species %in% multi_in_raw) / sum(sc$N == 1)
    t0 <- Sys.time(); r <- run_one(s$id, sp_map, sprintf("%s_s%d", arm, seed))
    el <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    cat(sprintf("N=%d seed=%d %-13s n=%5d n_spp=%5d singleton=%.1f%% (%.0f%% of them multi-seq in raw) | cut=%.3f F=%.4f F@1=%.4f edges=%d (%.0fs)\n",
                N_CAP, seed, arm, nrow(s), nrow(sc), 100 * mean(sc$N == 1), fake, r$cut, r$f, r$f_at_1, r$n_edges, el))
    out[[length(out) + 1]] <- data.table(N_cap = N_CAP, seed = seed, arm = arm, n = nrow(s), n_spp = nrow(sc),
      pct_singleton = round(100 * mean(sc$N == 1), 1), pct_singleton_multi_in_raw = round(fake, 1),
      cutoff = r$cut, f_measure = r$f, f_at_1 = r$f_at_1, n_edges = r$n_edges, elapsed_sec = el)
    fwrite(rbindlist(out), out_file, sep = "\t")
  }
  res <- rbindlist(out)
  cat("\n=== SUMMARY ===\n"); print(res[, .(cutoff_mean = mean(cutoff), cutoff_min = min(cutoff), cutoff_max = max(cutoff),
      pct_singleton = mean(pct_singleton), F = mean(f_measure)), by = arm])
}
