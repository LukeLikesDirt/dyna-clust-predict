#!/usr/bin/env Rscript
# diagnose_singletons.R — Quantify how singleton classes (species represented by
# a single sequence) drive species-level cut-offs toward 1.0 and distort the
# F-measure.
#
# Background: predict.R scores each class with a Dice term,
#   dice(c) = max_r 2*count(r,c) / (class_size(c) + cluster_size(r))
# weighted by class size. After dereplicate_lca.R every sequence in the dataset
# is unique, so at threshold 1.0 every sequence is its own cluster. A class of
# size 1 then scores 2*1/(1+1) = 1.0 -- a free perfect score that no lower
# threshold can beat, so any dataset made mostly of singleton species has its
# optimum pinned at 1.0 by construction.
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(data.table)
})

source("R/utils.R")

region     <- "full_ITS"
cutoffs_in <- file.path("data", region, "eukaryome.cutoffs.json.txt")
class_in   <- file.path("data", region, "eukaryome_ITS.classification")

cat("Region:", region, "\n\n")

cut_tbl <- read_tsv(cutoffs_in, show_col_types = FALSE) %>%
  rename(cutoff = `cut-off`, confidence = confidence,
         n_seq = `sequence number`, n_grp = `group number`)

cls <- fread(class_in, data.table = FALSE)

# ── Theoretical F-measure at threshold 1.0 ──────────────────────────────────
# Every sequence unique => every cluster is a singleton at t = 1.0.
#   class of size k  ->  best dice = 2/(k+1)
#   F(1.0) = sum_c [ k_c * 2/(k_c+1) ] / total_n
f_at_one <- function(k) sum(k * (2 / (k + 1))) / sum(k)

# ── Per-dataset species composition, for species-rank rows ──────────────────
species_rows <- cut_tbl %>% filter(rank == "species")

parent_col <- c(genus = "genus", family = "family", order = "order",
                class = "class", phylum = "phylum", kingdom = "kingdom")

setDT(cls)
cls_ok <- cls[is_identified(species)]

# One pass per parent rank: species counts within each parent taxon, then
# collapse to per-(parent taxon) composition stats. Avoids rescanning the
# full classification table once per dataset.
diag <- rbindlist(lapply(unique(species_rows$higher_rank), function(hr) {
  col <- parent_col[[hr]]
  if (is.null(col) || !col %in% names(cls_ok)) return(NULL)

  counts <- cls_ok[is_identified(get(col)),
                   .(k = .N), by = c(col, "species")]
  setnames(counts, col, "parent")

  stats <- counts[, .(
    n_species      = .N,
    n_seq_cls      = sum(k),
    singleton_frac = mean(k == 1),
    seq_per_sp     = mean(k),
    f_theory_at_1  = round(f_at_one(k), 4)
  ), by = parent]

  rows <- as.data.table(species_rows)[higher_rank == hr]
  merge(rows[, .(rank, higher_rank, dataset, cutoff, confidence)],
        stats, by.x = "dataset", by.y = "parent")
}), fill = TRUE)

cat("Species-rank datasets diagnosed:", nrow(diag), "\n\n")

# ── Does singleton fraction predict the cut-off? ────────────────────────────
cat("=== Cut-off vs singleton fraction ===\n")
diag %>%
  mutate(bin = cut(singleton_frac,
                   breaks = c(-0.01, 0.25, 0.5, 0.75, 0.95, 1.0),
                   labels = c("<25%", "25-50%", "50-75%", "75-95%", ">95%"))) %>%
  group_by(bin) %>%
  summarise(
    n              = n(),
    median_cutoff  = round(median(cutoff), 4),
    pct_cutoff_1   = round(100 * mean(cutoff >= 1), 1),
    pct_cutoff_999 = round(100 * mean(cutoff >= 0.999), 1),
    median_conf    = round(median(confidence), 3)
  ) %>% print(n = Inf)

cat("\nSpearman correlation singleton_frac vs cutoff: ",
    round(cor(diag$singleton_frac, diag$cutoff, method = "spearman"), 3), "\n")
cat("Spearman correlation singleton_frac vs confidence: ",
    round(cor(diag$singleton_frac, diag$confidence, method = "spearman"), 3), "\n")

# ── How well does the theoretical F(1.0) explain reported confidence? ───────
pinned <- diag %>% filter(cutoff >= 1)
cat("\n=== Datasets with cut-off pinned at 1.0 ===\n")
cat("n =", nrow(pinned), "\n")
cat("median singleton fraction:", round(median(pinned$singleton_frac), 3), "\n")
cat("median reported confidence:", round(median(pinned$confidence), 3), "\n")
cat("median THEORETICAL F at t=1.0 (from class sizes alone):",
    round(median(pinned$f_theory_at_1), 3), "\n")
cat("\nAgreement between reported confidence and theoretical F(1.0):\n")
cat("  mean |difference| =", round(mean(abs(pinned$confidence - pinned$f_theory_at_1)), 4), "\n")

cat("\nExamples (pinned at 1.0, sorted by singleton fraction):\n")
pinned %>%
  arrange(desc(singleton_frac)) %>%
  select(dataset, higher_rank, n_species, singleton_frac, confidence, f_theory_at_1) %>%
  head(12) %>% as.data.frame() %>% print(row.names = FALSE)

# ── How much F-measure weight do singletons carry overall? ──────────────────
cat("\n=== Share of F-measure weight held by singleton species ===\n")
diag %>%
  summarise(
    median_weight_share = round(median(singleton_frac * 1 / seq_per_sp), 3),
    mean_weight_share   = round(mean(singleton_frac * 1 / seq_per_sp), 3)
  ) %>% print()

dir.create("tests/CBSITS_eval/results", showWarnings = FALSE, recursive = TRUE)
write_tsv(diag, "tests/CBSITS_eval/results/singleton_diagnosis.tsv")
cat("\nWrote tests/CBSITS_eval/results/singleton_diagnosis.tsv\n")
