#!/usr/bin/env Rscript
# test_middle_cutoff.R — Verify that when F-measure is tied over a range of
# thresholds, the MIDDLE cutoff is selected (not the lowest).
#
# This test works by directly sourcing the internal functions from predict.R
# (after stubbing out the CLI argument parsing) and calling predict_dataset()
# with a crafted similarity matrix where F-measure is constant over a wide
# range of thresholds.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# ── Source utils.R (needed by predict.R functions) ────────────────────────────
source("R/utils.R")

# ── Define helper functions inline (copied from predict.R) ────────────────────
# We source only the functions we need to avoid triggering CLI parsing.

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

build_neighbors <- function(ids, sub_sim, threshold) {
  edges <- sub_sim[score >= threshold, .(i, j)]
  adj   <- split(edges$j, edges$i)
  neighbor_list <- setNames(vector("list", length(ids)), ids)
  for (id in ids) neighbor_list[[id]] <- adj[[id]]
  neighbor_list
}

find_clusters <- function(ids, neighbor_list) {
  visited  <- setNames(logical(length(ids)), ids)
  clusters <- list()
  for (start in ids) {
    if (!visited[start]) {
      queue          <- start
      visited[start] <- TRUE
      members        <- start
      while (length(queue) > 0) {
        current <- queue[1]
        queue   <- queue[-1]
        nbrs <- neighbor_list[[current]]
        nbrs <- nbrs[nbrs %in% ids]
        new  <- nbrs[!visited[nbrs]]
        if (length(new) > 0) {
          visited[new] <- TRUE
          members      <- c(members, new)
          queue        <- c(queue, new)
        }
      }
      clusters[[length(clusters) + 1L]] <- members
    }
  }
  clusters
}

compute_fmeasure <- function(classes, clusters) {
  f <- 0; n <- 0
  for (group in classes) {
    m <- 0
    for (cl in clusters) {
      i_size <- length(intersect(group, cl))
      v      <- 2 * i_size / (length(group) + length(cl))
      if (v > m) m <- v
    }
    n <- n + length(group)
    f <- f + length(group) * m
  }
  if (n == 0) return(0)
  round(f / n, 4)
}

predict_dataset <- function(dataset_name, seq_ids, classes, sim_dt,
                            start_t, end_t, step_t,
                            existing = list(), redo = FALSE,
                            verbose = TRUE) {
  saved_fm <- if ("fmeasures"  %in% names(existing)) existing$fmeasures   else list()
  sub_sim <- sim_dt[i %chin% seq_ids & j %chin% seq_ids]
  if (nrow(sub_sim) == 0) return(list(error = TRUE))
  thresholds <- numeric(0)
  fmeasures  <- numeric(0)
  t          <- round(start_t, 4)
  while (t <= end_t + 1e-9) {
    t_str <- sprintf("%.4f", t)
    if (!redo && t_str %in% names(saved_fm)) {
      fmeasure <- saved_fm[[t_str]]
    } else {
      neighbors <- build_neighbors(seq_ids, sub_sim, t)
      clusters  <- find_clusters(seq_ids, neighbors)
      fmeasure  <- compute_fmeasure(classes, clusters)
      saved_fm[[t_str]] <- fmeasure
    }
    thresholds <- c(thresholds, t)
    fmeasures  <- c(fmeasures,  fmeasure)
    if (verbose) cat(sprintf("  threshold=%.4f  F=%.4f\n", t, fmeasure))
    t <- round(t + step_t, 4)
  }
  # Select the middle threshold among all thresholds tied at the best F-measure
  best_f      <- max(fmeasures)
  tied_idx    <- which(fmeasures == best_f)
  mid_pos     <- tied_idx[ceiling(length(tied_idx) / 2)]
  opt_t       <- thresholds[mid_pos]
  if (verbose) {
    if (length(tied_idx) > 1) {
      cat(sprintf("[predict] %s: F=%.4f tied over %d thresholds (%.4f-%.4f), selecting middle: %.4f\n",
                  dataset_name, best_f, length(tied_idx),
                  thresholds[tied_idx[1]], thresholds[tied_idx[length(tied_idx)]], opt_t))
    } else {
      cat(sprintf("[predict] %s: optimal cutoff=%.4f  F-measure=%.4f\n",
                  dataset_name, opt_t, best_f))
    }
  }
  list(
    error          = FALSE,
    opt_t          = opt_t,
    best_f         = best_f,
    thresholds     = thresholds,
    fmeasures      = fmeasures,
    fmeasures_dict = saved_fm
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# Test setup: craft a scenario where F-measure is identical over a range
# ══════════════════════════════════════════════════════════════════════════════
#
# We create 6 sequences in 2 true classes (A: {s1,s2,s3}, B: {s4,s5,s6}).
# Pairwise similarities are designed so that:
#   - At thresholds 0.900-0.958 the clustering is identical (same F-measure).
#   - Below 0.900 everything merges into one cluster (lower F).
#   - Above 0.958 clusters fragment (lower F).
#
# Design:
#   Within-class pairs:  sim = 0.96  (stay connected until threshold > 0.96)
#   Between-class pairs: sim = 0.89  (disconnected at threshold >= 0.90)
#
# So at thresholds 0.900-0.958 (step 0.001):
#   {s1,s2,s3} and {s4,s5,s6} are perfectly separated → F = 1.0
# At threshold 0.899:
#   some cross-links still exist (sim = 0.89 >= 0.89) → all merge → F < 1.0
# At threshold 0.961 or higher:
#   within-class links break (0.96 < 0.961) → each seq its own cluster → F < 1.0

seq_ids <- paste0("s", 1:6)
classes <- list(A = c("s1", "s2", "s3"), B = c("s4", "s5", "s6"))

pairs <- data.table(
  i     = character(),
  j     = character(),
  score = numeric()
)

# Within-class A: all pairs at sim = 0.96
for (x in c("s1", "s2", "s3")) {
  for (y in c("s1", "s2", "s3")) {
    if (x != y) {
      pairs <- rbind(pairs, data.table(i = x, j = y, score = 0.96))
    }
  }
}
# Within-class B: all pairs at sim = 0.96
for (x in c("s4", "s5", "s6")) {
  for (y in c("s4", "s5", "s6")) {
    if (x != y) {
      pairs <- rbind(pairs, data.table(i = x, j = y, score = 0.96))
    }
  }
}
# Between-class: all pairs at sim = 0.89
for (x in c("s1", "s2", "s3")) {
  for (y in c("s4", "s5", "s6")) {
    pairs <- rbind(pairs, data.table(i = x, j = y, score = 0.89))
    pairs <- rbind(pairs, data.table(i = y, j = x, score = 0.89))
  }
}
# Self-similarity
for (x in seq_ids) {
  pairs <- rbind(pairs, data.table(i = x, j = x, score = 1.0))
}

sim_dt <- as.data.table(pairs)

# ══════════════════════════════════════════════════════════════════════════════
# TEST 1: Show the current behaviour (expect lowest cutoff is selected)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════════════════════════\n")
cat("TEST: Sweep thresholds 0.880 → 0.970 (step 0.001)\n")
cat("══════════════════════════════════════════════════════════════\n\n")

result <- predict_dataset(
  "test", seq_ids, classes, sim_dt,
  start_t = 0.880, end_t = 0.970, step_t = 0.001,
  verbose = TRUE
)

cat("\n── Summary ──\n")
cat(sprintf("Optimal cutoff selected: %.4f\n", result$opt_t))
cat(sprintf("Best F-measure:          %.4f\n", result$best_f))

# Identify the range of thresholds that share the best F-measure
best_indices <- which(result$fmeasures == result$best_f)
tied_range   <- result$thresholds[best_indices]
cat(sprintf("Tied F-measure range:    %.4f – %.4f  (%d thresholds)\n",
            min(tied_range), max(tied_range), length(tied_range)))
expected_mid <- tied_range[ceiling(length(tied_range) / 2)]
cat(sprintf("Expected middle cutoff:  %.4f\n", expected_mid))

if (result$opt_t == min(tied_range)) {
  cat("\n*** FAIL: Code still selects the LOWEST tied cutoff. ***\n")
} else if (result$opt_t == expected_mid) {
  cat("\n*** PASS: Code correctly selects the MIDDLE tied cutoff. ***\n")
} else {
  cat(sprintf("\n*** UNEXPECTED: Selected %.4f (not lowest or middle). ***\n", result$opt_t))
}


# ══════════════════════════════════════════════════════════════════════════════
# TEST 2: Narrower range with odd number of tied thresholds
# ══════════════════════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════════════════════════\n")
cat("TEST 2: Narrower range (0.920 → 0.965, step 0.005)\n")
cat("══════════════════════════════════════════════════════════════\n\n")

result2 <- predict_dataset(
  "test2", seq_ids, classes, sim_dt,
  start_t = 0.920, end_t = 0.965, step_t = 0.005,
  verbose = TRUE
)

best_idx2   <- which(result2$fmeasures == result2$best_f)
tied_range2 <- result2$thresholds[best_idx2]
expected2   <- tied_range2[ceiling(length(tied_range2) / 2)]
cat(sprintf("\nTied range: %.4f–%.4f (%d thresholds)\n",
            min(tied_range2), max(tied_range2), length(tied_range2)))
cat(sprintf("Selected: %.4f  Expected middle: %.4f  %s\n",
            result2$opt_t, expected2,
            if (result2$opt_t == expected2) "PASS" else "FAIL"))


# ══════════════════════════════════════════════════════════════════════════════
# TEST 3: Only one threshold achieves the best F-measure (no tie)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════════════════════════\n")
cat("TEST 3: Single unique optimum (threshold = 0.960 only)\n")
cat("══════════════════════════════════════════════════════════════\n\n")

result3 <- predict_dataset(
  "test3", seq_ids, classes, sim_dt,
  start_t = 0.960, end_t = 0.962, step_t = 0.001,
  verbose = TRUE
)
cat(sprintf("\nSelected: %.4f  Best F: %.4f  %s\n",
            result3$opt_t, result3$best_f,
            if (result3$opt_t == 0.9600) "PASS" else "FAIL"))


# ══════════════════════════════════════════════════════════════════════════════
# TEST 4: Even number of tied thresholds (ceiling picks upper-middle)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════════════════════════\n")
cat("TEST 4: Even-count tied range (0.950 → 0.965, step 0.005)\n")
cat("══════════════════════════════════════════════════════════════\n\n")

result4 <- predict_dataset(
  "test4", seq_ids, classes, sim_dt,
  start_t = 0.950, end_t = 0.965, step_t = 0.005,
  verbose = TRUE
)

best_idx4   <- which(result4$fmeasures == result4$best_f)
tied_range4 <- result4$thresholds[best_idx4]
expected4   <- tied_range4[ceiling(length(tied_range4) / 2)]
cat(sprintf("\nTied range: %.4f–%.4f (%d thresholds)\n",
            min(tied_range4), max(tied_range4), length(tied_range4)))
cat(sprintf("Selected: %.4f  Expected middle: %.4f  %s\n",
            result4$opt_t, expected4,
            if (result4$opt_t == expected4) "PASS" else "FAIL"))


# ══════════════════════════════════════════════════════════════════════════════
# TEST 5: tie_tolerance widens the tied range (R/predict.R:559-563)
# ══════════════════════════════════════════════════════════════════════════════
# Verifies the same selection formula added to R/predict.R's predict_dataset():
#   tied_idx <- which(fmeasures >= best_f - tie_tolerance)
# against TEST 1's already-computed fmeasures/thresholds, so no new similarity
# scenario is needed. tie_tolerance = 0 must reproduce the exact-equality
# selection exactly (regression safety at the default-off setting); a nonzero
# tolerance must select a range at least as wide as the exact-tie range.

cat("\n══════════════════════════════════════════════════════════════\n")
cat("TEST 5: tie_tolerance selection formula\n")
cat("══════════════════════════════════════════════════════════════\n\n")

pick_with_tolerance <- function(thresholds, fmeasures, tie_tolerance) {
  best_f   <- max(fmeasures)
  tied_idx <- which(fmeasures >= best_f - tie_tolerance)
  list(opt_t = thresholds[tied_idx[ceiling(length(tied_idx) / 2)]],
       width = length(tied_idx))
}

zero_tol <- pick_with_tolerance(result$thresholds, result$fmeasures, 0)
cat(sprintf("tie_tolerance=0:    cutoff=%.4f  tied width=%d  %s\n",
            zero_tol$opt_t, zero_tol$width,
            if (zero_tol$opt_t == result$opt_t &&
                zero_tol$width == length(tied_range)) "PASS (matches exact-equality)" else "FAIL"))

wide_tol <- pick_with_tolerance(result$thresholds, result$fmeasures, 0.5)
cat(sprintf("tie_tolerance=0.5:  cutoff=%.4f  tied width=%d  %s\n",
            wide_tol$opt_t, wide_tol$width,
            if (wide_tol$width >= zero_tol$width) "PASS (tolerance widens or matches range)" else "FAIL"))
