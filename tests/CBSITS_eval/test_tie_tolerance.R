#!/usr/bin/env Rscript
# test_tie_tolerance.R — Phase 1d: how much does the cut-off selection rule
# change if ties are widened from exact-equality to a tolerance band?
#
# predict.R currently selects: tied_idx <- which(fmeasures == best_f), where
# fmeasures is rounded to 4 dp (R/predict.R:520,560). This re-derives the
# optimum from the SAME cached per-threshold F traces already stored in
# <prefix>.predicted, under alternative rules:
#   - current:  round to 4dp, exact equality
#   - round3:   round to 3dp, exact equality
#   - tol=t:    4dp values, tied_idx <- which(fmeasures >= best_f - t)
# No vsearch, no recompute of similarity -- this is free.
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
})

regions <- c("full_ITS", "ITS2")
tolerances <- c(0.0005, 0.001, 0.002)

pick_current <- function(thresholds, fm) {
  best <- max(fm)
  tied <- which(fm == best)
  list(cut = thresholds[tied[ceiling(length(tied) / 2)]], f = best, width = length(tied))
}

pick_round3 <- function(thresholds, fm) {
  fm3 <- round(fm, 3)
  best <- max(fm3)
  tied <- which(fm3 == best)
  list(cut = thresholds[tied[ceiling(length(tied) / 2)]], f = max(fm[tied]), width = length(tied))
}

pick_tol <- function(thresholds, fm, tol) {
  best <- max(fm)
  tied <- which(fm >= best - tol)
  # tied need not be contiguous; use full range for "width" reporting
  list(cut = thresholds[tied[ceiling(length(tied) / 2)]], f = best, width = length(tied))
}

process_region <- function(region) {
  path <- file.path("data", region, "eukaryome.predicted")
  cat("\n==============================\nRegion:", region, "\n==============================\n")
  cat("Loading:", path, "...\n")
  d <- fromJSON(path, simplifyVector = FALSE)

  rows <- list()
  for (rank in names(d)) {
    datasets <- d[[rank]]
    for (dname in names(datasets)) {
      ds <- datasets[[dname]]
      fm_list <- ds[["fmeasures"]]
      if (is.null(fm_list) || length(fm_list) == 0) next

      t_str <- names(fm_list)
      thresholds <- as.numeric(t_str)
      fm <- as.numeric(unlist(fm_list))
      ord <- order(thresholds)
      thresholds <- thresholds[ord]; fm <- fm[ord]

      cur <- pick_current(thresholds, fm)
      r3  <- pick_round3(thresholds, fm)
      t05 <- pick_tol(thresholds, fm, 0.0005)
      t10 <- pick_tol(thresholds, fm, 0.001)
      t20 <- pick_tol(thresholds, fm, 0.002)

      rows[[length(rows) + 1]] <- data.table(
        region = region, rank = rank, higher_rank = ds[["higher_rank"]] %||% NA,
        dataset = dname,
        n_seq = ds[["sequence number"]] %||% NA, n_grp = ds[["group number"]] %||% NA,
        cut_current = cur$cut, f_current = cur$f, width_current = cur$width,
        cut_round3  = r3$cut,  f_round3  = r3$f,  width_round3  = r3$width,
        cut_tol0005 = t05$cut, width_tol0005 = t05$width,
        cut_tol001  = t10$cut, width_tol001  = t10$width,
        cut_tol002  = t20$cut, width_tol002  = t20$width
      )
    }
  }
  rbindlist(rows)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

all_res <- rbindlist(lapply(regions, process_region))

cat("\n\n=== SUMMARY: how many cut-offs move, by rule ===\n")
for (col in c("round3", "tol0005", "tol001", "tol002")) {
  cutcol <- paste0("cut_", col)
  d <- copy(all_res)
  d[, delta := abs(get(cutcol) - cut_current)]
  cat(sprintf("\n--- %s ---\n", col))
  print(d[, .(
    n = .N,
    n_changed = sum(delta > 1e-9),
    pct_changed = round(100 * mean(delta > 1e-9), 1),
    median_abs_delta = round(median(delta), 4),
    mean_abs_delta = round(mean(delta), 4),
    max_abs_delta = round(max(delta), 4)
  ), by = .(region, rank)][order(region, factor(rank, levels = c("species","genus","family","order","class","phylum")))])
}

cat("\n\n=== Tie-width distribution (current rule, 4dp exact) ===\n")
print(all_res[, .(
  n = .N,
  median_width = as.numeric(median(width_current)),
  p90_width = as.numeric(quantile(width_current, 0.9)),
  max_width = as.numeric(max(width_current))
), by = .(region, rank)][order(region, factor(rank, levels = c("species","genus","family","order","class","phylum")))])

cat("\n\n=== Pinned-at-1.0 rows: do they move under any rule? ===\n")
pinned <- all_res[rank == "species" & cut_current >= 1]
cat("n pinned at 1.0 (species rank):", nrow(pinned), "\n")
for (col in c("round3", "tol0005", "tol001", "tol002")) {
  cutcol <- paste0("cut_", col)
  n_moved <- pinned[get(cutcol) < 1, .N]
  cat(sprintf("  moved off 1.0 under %-8s: %d / %d\n", col, n_moved, nrow(pinned)))
}

dir.create("tests/CBSITS_eval/results", showWarnings = FALSE, recursive = TRUE)
fwrite(all_res, "tests/CBSITS_eval/results/tie_tolerance_test.tsv", sep = "\t")
cat("\nWrote tests/CBSITS_eval/results/tie_tolerance_test.tsv (", nrow(all_res), "rows )\n")
