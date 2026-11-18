#!/usr/bin/env Rscript
# selection_test.R — Compare structured vs random sequence selection
#
# Creates two subsets from a larger test dataset:
#   1. Structured: hierarchical budget allocation (subset.R logic)
#   2. Random: simple random.sample (dnabarcoder logic)
#
# Then runs predict.R on both and compares cutoffs/confidence.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

source("R/utils.R")

# ── Parameters ────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

input_fasta   <- args[1]  # test_selection.fasta
input_class   <- args[2]  # test_selection.classification
output_dir    <- args[3]  # tests/comparison_test/data
max_sequences <- as.integer(args[4])  # e.g. 500

rank_hierarchy <- c("kingdom", "phylum", "class", "order", "family", "genus")

cat("=== SELECTION TEST ===\n")
cat(sprintf("Input FASTA:   %s\n", input_fasta))
cat(sprintf("Classification: %s\n", input_class))
cat(sprintf("Max sequences:  %d\n", max_sequences))
cat(sprintf("Output dir:     %s\n", output_dir))
cat("\n")

# ── Read classification ──────────────────────────────────────────────────────
cls <- fread(input_class, sep = "\t", header = TRUE, quote = "",
             fill = TRUE, na.strings = "", check.names = FALSE,
             data.table = FALSE)

# ── Functions (from subset.R) ─────────────────────────────────────────────────

balanced_downsample <- function(df, target_rank, n) {
  df %>%
    ungroup() %>%
    group_by(!!sym(target_rank)) %>%
    mutate(.draw_order = row_number()) %>%
    ungroup() %>%
    arrange(.draw_order, !!sym(target_rank)) %>%
    slice_head(n = n) %>%
    select(-.draw_order)
}

allocate_budget <- function(available, budget) {
  n <- length(available)
  if (n == 0 || budget == 0) return(setNames(integer(n), names(available)))
  if (sum(available) <= budget) return(available)

  allocation <- setNames(integer(n), names(available))
  remaining  <- budget
  competing  <- rep(TRUE, n)

  while (remaining > 0 && any(competing)) {
    n_competing <- sum(competing)
    per_group   <- remaining / n_competing
    underfull   <- competing & (available < per_group)

    if (!any(underfull)) {
      base_share              <- floor(remaining / n_competing)
      allocation[competing]  <- base_share
      leftover               <- remaining - base_share * n_competing
      if (leftover > 0) {
        top_up_idx              <- which(competing)[seq_len(leftover)]
        allocation[top_up_idx] <- allocation[top_up_idx] + 1L
      }
      break
    }
    for (i in which(underfull)) {
      allocation[i] <- available[i]
      remaining     <- remaining - available[i]
      competing[i]  <- FALSE
    }
  }
  allocation
}

hierarchical_sample <- function(df, grouping_ranks, target_rank, budget) {
  if (nrow(df) == 0 || budget == 0) return(NULL)
  if (length(grouping_ranks) == 0) {
    return(balanced_downsample(df, target_rank, min(budget, nrow(df))))
  }

  current_rank    <- grouping_ranks[1]
  remaining_ranks <- grouping_ranks[-1]
  df              <- df %>% filter(is_identified(!!sym(current_rank)))
  if (nrow(df) == 0) return(NULL)

  groups      <- split(df, df[[current_rank]])
  available   <- sapply(groups, nrow)
  allocations <- allocate_budget(available, budget)

  results <- list()
  for (grp in names(groups)) {
    if (allocations[[grp]] == 0) next
    res <- hierarchical_sample(groups[[grp]], remaining_ranks, target_rank, allocations[[grp]])
    if (!is.null(res) && nrow(res) > 0) results[[grp]] <- res
  }
  if (length(results) == 0) return(NULL)
  bind_rows(results)
}

global_prediction_filter <- function(df, target_rank, hierarchy = rank_hierarchy,
                                     max_sequences = 25000, min_sequences = 30) {
  target_idx     <- which(hierarchy == target_rank)
  grouping_ranks <- hierarchy[seq_len(target_idx - 1)]

  df <- df %>% filter(is_identified(!!sym(target_rank)))
  if (nrow(df) < min_sequences) return(NULL)

  if (length(grouping_ranks) == 0) {
    result <- balanced_downsample(df, target_rank, min(max_sequences, nrow(df)))
    if (nrow(result) < min_sequences) return(NULL)
    return(result)
  }

  result <- hierarchical_sample(df, grouping_ranks, target_rank, max_sequences)
  if (is.null(result) || nrow(result) < min_sequences) return(NULL)
  result
}

# ── For each rank: produce structured and random ID lists ─────────────────────

for (rank in rank_hierarchy) {
  cat(sprintf("── Rank: %s ──\n", rank))

  # Filter to identified sequences at this rank
  valid <- cls[is_identified(cls[[rank]]), ]
  cat(sprintf("  Identified sequences: %d\n", nrow(valid)))

  if (nrow(valid) < 30) {
    cat("  Too few sequences, skipping.\n\n")
    next
  }

  # ── Structured selection ──────────────────────────────────────────────────
  set.seed(42)
  structured <- global_prediction_filter(valid, rank,
                                         hierarchy = rank_hierarchy,
                                         max_sequences = max_sequences,
                                         min_sequences = 30)
  if (!is.null(structured)) {
    structured_ids <- structured$id
    cat(sprintf("  Structured: %d sequences selected\n", length(structured_ids)))

    # Report group distribution
    grp_counts <- table(structured[[rank]])
    cat(sprintf("    Groups: %d, max proportion: %.3f\n",
                length(grp_counts), max(grp_counts) / sum(grp_counts)))
  } else {
    cat("  Structured: filtered out (not enough data)\n")
    structured_ids <- NULL
  }

  # ── Random selection (dnabarcoder style) ──────────────────────────────────
  set.seed(42)
  if (nrow(valid) > max_sequences) {
    random_ids <- sample(valid$id, max_sequences)
  } else {
    random_ids <- valid$id
  }
  cat(sprintf("  Random:     %d sequences selected\n", length(random_ids)))

  random_valid <- valid[valid$id %in% random_ids, ]
  grp_counts_r <- table(random_valid[[rank]])
  cat(sprintf("    Groups: %d, max proportion: %.3f\n",
              length(grp_counts_r), max(grp_counts_r) / sum(grp_counts_r)))

  # ── Save ID lists ─────────────────────────────────────────────────────────
  if (!is.null(structured_ids)) {
    writeLines(structured_ids,
               file.path(output_dir, sprintf("selection_structured_%s_ids.txt", rank)))
  }
  writeLines(random_ids,
             file.path(output_dir, sprintf("selection_random_%s_ids.txt", rank)))

  cat("\n")
}

cat("ID lists saved. Now run predict.R on each subset.\n")
