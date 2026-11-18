#!/usr/bin/env Rscript
# test_singleton_constraints.R — Test two ideas that act directly on the
# sequence POOL before the F-measure sweep runs, rather than patching the
# result afterward (0.995 cap) or rebalancing clade proportions (cross-rank
# max_proportion). Both target the actual mechanism identified when
# explaining current singleton handling: a species with exactly one sequence
# scores a free Dice=1.0 at threshold=1.0 (2*1/(1+1)), and neither existing
# guard (--min_multiseq_groups eligibility gate in R/predict.R, the
# consolidation override guard in R/consolidate_cutoffs.R) stops singleton
# species from contributing that free score to the weighted F total once a
# dataset clears the eligibility bar -- they only gate whether/how a result
# gets reported, not what goes into computing it.
#
#   (b) max 50% singleton-species constraint: if singleton species (1
#       sequence) outnumber multi-sequence species (>=2 sequences), randomly
#       drop excess singleton species until they're at most 50% of the
#       species COUNT (not sequence count) in the pool.
#   (c) min 2 sequences per group: drop every singleton species from the
#       pool entirely before the sweep ever runs.
#
# Both are applied on top of the existing species-dominance cap
# (R/subset.R's cap_dominant_group, max_proportion=0.5) so the comparison
# isolates the new constraint's effect rather than also changing existing
# behaviour. (a) current/baseline is reused from
# tests/CBSITS_eval/results/cross_rank_balance.tsv for the 11 groups already
# tested there (identical pool-construction + N_CAP=5000 + seed=42), and
# computed fresh for 3 new genus-level mega-genus groups (Aspergillus,
# Fusarium, Colletotrichum) added here to check the effect on genuinely
# species-complex-driven (not singleton-dominated) pinned groups.
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
tmp_dir   <- "tmp/singleton_constraints_test"
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

start_t <- 0.90; end_t <- 1.0; step_t <- 0.001
N_CAP <- 5000L
MAX_PROPORTION <- 0.5
MAX_SINGLETON_FRAC <- 0.5
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
cat("Testing", nrow(targets), "(parent_rank, dataset) pairs:\n")
print(as.data.frame(targets))

cls <- fread(class_all, sep = "\t", header = TRUE, colClasses = "character")
cls <- cls[is_identified(species)]

# Reuse baseline "current" results already computed under identical
# pool-construction + N_CAP + seed for the 11 groups shared with the
# cross-rank balance test.
baseline <- fread("tests/CBSITS_eval/results/cross_rank_balance.tsv", sep = "\t", header = TRUE)
baseline <- baseline[, .(parent_rank, dataset, n_current, n_spp_current, cut_current, f_current)]

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

# Randomly drop excess singleton (1-sequence) species until they're at most
# max_singleton_frac of the species COUNT (not sequence count).
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

# Drop every singleton species from the pool entirely.
drop_singletons <- function(df) {
  spp_counts <- df[, .N, by = species]
  multi_spp <- spp_counts[N >= 2, species]
  df[species %in% multi_spp]
}

# Downsample to at most n_cap sequences by dropping WHOLE species at random,
# never truncating a kept species' own sequences. A plain row-level
# slice_sample() (what R/subset.R's Filter 4 actually does in production)
# silently undoes cap_singleton_fraction()/drop_singletons() above: sampling
# individual sequences uniformly fragments most multi-sequence species down
# to exactly one surviving sequence once the sampling fraction is steep
# (e.g. 5000/60000 ~ 8%), turning them back into apparent singletons in the
# final pool -- confirmed by a first run of this script where the min2seq
# arm (which should be 0% singleton by construction) came out 84.9%
# singleton for kingdom Fungi after a naive row-level downsample. Species-
# level downsampling is the only way to honestly preserve either
# constraint through to the pool vsearch actually sees.
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

# ── F-measure sweep (same union-find pattern as prior CBSITS_eval scripts) ──

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
  n_single <- sum(spp_counts$N == 1)
  100 * n_single / nrow(spp_counts)
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

  # Baseline pool (species-dominance-capped, pre-downsample) -- needed as the
  # starting point for both new constraints, applied ON TOP of it.
  base_pool <- cap_dominant_group(pool, "species", MAX_PROPORTION)

  # (b) max 50% singleton-species fraction, then downsample by whole species
  #     (never fragment a kept species below its actual sequence count)
  b_pool <- cap_singleton_fraction(base_pool, MAX_SINGLETON_FRAC)
  b_pool <- downsample_by_species(b_pool, N_CAP)

  # (c) drop all singletons, then downsample by whole species
  c_pool <- drop_singletons(base_pool)
  c_pool <- downsample_by_species(c_pool, N_CAP)

  cat(sprintf("  max50pct_singleton: n=%d, %.1f%% singleton spp, %d spp\n",
              nrow(b_pool), singleton_pct(b_pool), length(unique(b_pool$species))))
  cat(sprintf("  min2seq           : n=%d, %.1f%% singleton spp, %d spp\n",
              nrow(c_pool), singleton_pct(c_pool), length(unique(c_pool$species))))

  n_threads <- if (max(nrow(b_pool), nrow(c_pool)) > 500) 10L else 4L
  safe <- gsub("[^A-Za-z0-9_]", "_", paste(pr, ds, sep = "_"))

  t0 <- Sys.time()
  r_b <- if (nrow(b_pool) >= 30) run_vsearch_f(b_pool$id, sp_map, paste0(safe, "_b"), n_threads) else list(cut = NA, f = NA)
  r_c <- if (nrow(c_pool) >= 30) run_vsearch_f(c_pool$id, sp_map, paste0(safe, "_c"), n_threads) else list(cut = NA, f = NA)

  # (a) baseline: reuse from cross_rank_balance.tsv if present, else compute fresh
  bl <- baseline[parent_rank == pr & dataset == ds]
  if (nrow(bl) == 1) {
    r_a <- list(cut = bl$cut_current, f = bl$f_current)
    n_a <- bl$n_current; n_spp_a <- bl$n_spp_current
  } else {
    a_pool <- base_pool
    if (nrow(a_pool) > N_CAP) a_pool <- a_pool[sample(.N, N_CAP)]
    r_a <- run_vsearch_f(a_pool$id, sp_map, paste0(safe, "_a"), n_threads)
    n_a <- nrow(a_pool); n_spp_a <- length(unique(a_pool$species))
  }
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)

  cat(sprintf("  [%d/%d] (%.1fs) current: cut=%s F=%s | max50pct_singleton: cut=%s F=%s | min2seq: cut=%s F=%s\n",
              i, nrow(targets), elapsed,
              ifelse(is.na(r_a$cut), "NA", sprintf("%.3f", r_a$cut)), ifelse(is.na(r_a$f), "NA", sprintf("%.4f", r_a$f)),
              ifelse(is.na(r_b$cut), "NA", sprintf("%.3f", r_b$cut)), ifelse(is.na(r_b$f), "NA", sprintf("%.4f", r_b$f)),
              ifelse(is.na(r_c$cut), "NA", sprintf("%.3f", r_c$cut)), ifelse(is.na(r_c$f), "NA", sprintf("%.4f", r_c$f))))

  out[[length(out) + 1]] <- data.table(
    parent_rank = pr, dataset = ds,
    n_current = n_a, n_spp_current = n_spp_a, cut_current = r_a$cut, f_current = r_a$f,
    n_max50 = nrow(b_pool), n_spp_max50 = length(unique(b_pool$species)),
    pct_singleton_max50 = round(singleton_pct(b_pool), 1),
    cut_max50 = r_b$cut, f_max50 = r_b$f,
    n_min2seq = nrow(c_pool), n_spp_min2seq = length(unique(c_pool$species)),
    cut_min2seq = r_c$cut, f_min2seq = r_c$f,
    elapsed_sec = elapsed
  )
  fwrite(rbindlist(out), "tests/CBSITS_eval/results/singleton_constraints.tsv", sep = "\t")
}

res <- rbindlist(out)
cat("\n=== SUMMARY ===\n")
print(as.data.frame(res), row.names = FALSE)
cat("\nWrote tests/CBSITS_eval/results/singleton_constraints.tsv\n")
