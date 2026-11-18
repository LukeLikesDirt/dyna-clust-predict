#!/usr/bin/env Rscript
# selection_test_4way.R — Compare four sequence selection approaches
#
# 1. random:          Pure random sample (dnabarcoder style)
# 2. random_capped:   Random sample but dominant group capped at ≤50%
# 3. structured:      Hierarchical budget allocation (subset.R logic)
# 4. structured_min5: Hierarchical but each subgroup must have ≥5 seqs
#
# Usage:
#   Rscript selection_test_4way.R <fasta> <classification> <output_dir> <max_seqs>

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

source("R/utils.R")

args <- commandArgs(trailingOnly = TRUE)
input_fasta   <- args[1]
input_class   <- args[2]
output_dir    <- args[3]
max_sequences <- as.integer(args[4])
max_prop      <- 0.5

rank_hierarchy <- c("kingdom", "phylum", "class", "order", "family", "genus")

cat("=== 4-WAY SELECTION TEST ===\n")
cat(sprintf("Max sequences: %d   Max proportion: %.2f\n\n", max_sequences, max_prop))

cls <- fread(input_class, sep = "\t", header = TRUE, quote = "",
             fill = TRUE, na.strings = "", check.names = FALSE,
             data.table = FALSE)

# ── Helpers from subset.R ─────────────────────────────────────────────────────

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
      base_share             <- floor(remaining / n_competing)
      allocation[competing] <- base_share
      leftover              <- remaining - base_share * n_competing
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

hierarchical_sample <- function(df, grouping_ranks, target_rank, budget,
                                min_per_group = 0L) {
  if (nrow(df) == 0 || budget == 0) return(NULL)
  if (length(grouping_ranks) == 0) {
    if (min_per_group > 0L) {
      # Remove groups with fewer than min_per_group sequences in source
      grp_sizes <- df %>% count(!!sym(target_rank)) %>% arrange(n)
      keep <- grp_sizes %>% filter(n >= min_per_group) %>% pull(!!sym(target_rank))
      df <- df %>% filter(!!sym(target_rank) %in% keep)
      if (nrow(df) == 0) return(NULL)
      n_groups <- length(keep)
      # Also ensure budget can give each remaining group >= min_per_group
      while (n_groups > 0 && budget < n_groups * min_per_group) {
        # Drop the smallest group to make room
        grp_sizes <- grp_sizes %>% filter(!!sym(target_rank) %in% keep)
        keep <- keep[-1]  # drop smallest
        n_groups <- length(keep)
        df <- df %>% filter(!!sym(target_rank) %in% keep)
      }
      if (nrow(df) == 0 || n_groups == 0) return(NULL)
    }
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
    res <- hierarchical_sample(groups[[grp]], remaining_ranks, target_rank,
                               allocations[[grp]], min_per_group)
    if (!is.null(res) && nrow(res) > 0) results[[grp]] <- res
  }
  if (length(results) == 0) return(NULL)
  bind_rows(results)
}

# ── Selection approach 1: Pure random ─────────────────────────────────────────

select_random <- function(valid, rank, max_n) {
  if (nrow(valid) <= max_n) return(valid$id)
  sample(valid$id, max_n)
}

# ── Selection approach 2: Random + cap dominant group at ≤50% ─────────────────

select_random_capped <- function(valid, rank, max_n, max_prop = 0.5) {
  grp_counts <- valid %>% count(!!sym(rank), sort = TRUE)
  dom_name   <- grp_counts[[rank]][1]
  dom_ids    <- valid %>% filter(!!sym(rank) == dom_name) %>% pull(id)
  rest_ids   <- valid %>% filter(!!sym(rank) != dom_name) %>% pull(id)

  if (length(rest_ids) == 0) return(sample(dom_ids, min(max_n, length(dom_ids))))

  dom_frac <- length(dom_ids) / nrow(valid)
  if (dom_frac <= max_prop) {
    # No cap needed
    if (nrow(valid) <= max_n) return(valid$id)
    return(sample(valid$id, max_n))
  }

  # Determine rest budget and dominant cap
  rest_avail <- length(rest_ids)
  if (rest_avail >= floor(max_n * (1 - max_prop))) {
    # Enough rest to fill (1-max_prop) of max_n
    rest_n <- floor(max_n * (1 - max_prop))
    dom_n  <- max_n - rest_n
  } else {
    # Not enough rest: take all rest, cap dominant at rest count
    rest_n <- rest_avail
    dom_n  <- rest_n  # so dominant = rest → 50/50
  }
  dom_n  <- min(dom_n, length(dom_ids))
  rest_n <- min(rest_n, rest_avail)

  c(sample(dom_ids, dom_n), sample(rest_ids, rest_n))
}

# ── Selection approach 3: Structured hierarchical ─────────────────────────────

select_structured <- function(valid, rank, max_n, hierarchy) {
  target_idx     <- which(hierarchy == rank)
  grouping_ranks <- hierarchy[seq_len(target_idx - 1)]

  df <- valid %>% filter(is_identified(!!sym(rank)))
  if (nrow(df) < 30) return(NULL)

  if (length(grouping_ranks) == 0) {
    result <- balanced_downsample(df, rank, min(max_n, nrow(df)))
    return(result$id)
  }

  result <- hierarchical_sample(df, grouping_ranks, rank, max_n,
                                min_per_group = 0L)
  if (is.null(result) || nrow(result) < 30) return(NULL)
  result$id
}

# ── Selection approach 4: Structured + min 5 per group ────────────────────────

select_structured_min5 <- function(valid, rank, max_n, hierarchy) {
  target_idx     <- which(hierarchy == rank)
  grouping_ranks <- hierarchy[seq_len(target_idx - 1)]

  df <- valid %>% filter(is_identified(!!sym(rank)))
  if (nrow(df) < 30) return(NULL)

  if (length(grouping_ranks) == 0) {
    # Remove rare groups first
    grp_sizes <- df %>% count(!!sym(rank))
    keep <- grp_sizes %>% filter(n >= 5) %>% pull(!!sym(rank))
    df <- df %>% filter(!!sym(rank) %in% keep)
    if (nrow(df) < 30) return(NULL)
    result <- balanced_downsample(df, rank, min(max_n, nrow(df)))
    return(result$id)
  }

  result <- hierarchical_sample(df, grouping_ranks, rank, max_n,
                                min_per_group = 5L)
  if (is.null(result) || nrow(result) < 30) return(NULL)
  result$id
}

# ── Run all four approaches for each rank ────────────────────────────────────

approaches <- c("random", "random_capped", "structured", "structured_min5")

for (rank in rank_hierarchy) {
  cat(sprintf("── Rank: %s ──\n", rank))

  valid <- cls[is_identified(cls[[rank]]), ]
  cat(sprintf("  Total identified: %d seqs, %d groups\n",
              nrow(valid), length(unique(valid[[rank]]))))

  if (nrow(valid) < 30) {
    cat("  Skipping (too few).\n\n")
    next
  }

  for (approach in approaches) {
    set.seed(42)

    ids <- switch(approach,
      "random"           = select_random(valid, rank, max_sequences),
      "random_capped"    = select_random_capped(valid, rank, max_sequences, max_prop),
      "structured"       = select_structured(valid, rank, max_sequences, rank_hierarchy),
      "structured_min5"  = select_structured_min5(valid, rank, max_sequences, rank_hierarchy)
    )

    if (is.null(ids) || length(ids) < 30) {
      cat(sprintf("  %-20s: FILTERED OUT\n", approach))
      next
    }

    sub <- valid[valid$id %in% ids, ]
    grp_counts <- table(sub[[rank]])
    n_groups   <- length(grp_counts)
    max_p      <- round(max(grp_counts) / sum(grp_counts), 3)
    min_grp_sz <- min(grp_counts)

    cat(sprintf("  %-20s: %4d seqs  %4d groups  max_prop=%.3f  min_grp_size=%d\n",
                approach, length(ids), n_groups, max_p, min_grp_sz))

    writeLines(ids,
               file.path(output_dir,
                         sprintf("sel4_%s_%s_ids.txt", approach, rank)))
  }
  cat("\n")
}

cat("All ID lists saved.\n")
