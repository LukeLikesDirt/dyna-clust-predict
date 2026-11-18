#!/usr/bin/env Rscript
# compare_cutoffs.R — Compare CBSITS (curated, strain/ex-type-vouchered) ITS2
# cutoff predictions against the existing EUKARYOME ITS2 cutoffs, for every
# (rank, higher_rank, dataset) triple predicted in both.
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

euk_path <- "data/ITS2/eukaryome.cutoffs.json.txt"
cbs_path <- "data/CBSITS_ITS2/cbsits.cutoffs.json.txt"

euk <- read_tsv(euk_path, show_col_types = FALSE) %>%
  rename(
    euk_cutoff = `cut-off`, euk_confidence = confidence,
    euk_n_seq = `sequence number`, euk_n_group = `group number`,
    euk_max_prop = `max proportion`
  )

cbs <- read_tsv(cbs_path, show_col_types = FALSE) %>%
  rename(
    cbs_cutoff = `cut-off`, cbs_confidence = confidence,
    cbs_n_seq = `sequence number`, cbs_n_group = `group number`,
    cbs_max_prop = `max proportion`
  )

joined <- inner_join(cbs, euk, by = c("rank", "higher_rank", "dataset")) %>%
  mutate(
    diff_cutoff     = round(cbs_cutoff - euk_cutoff, 4),
    diff_confidence = round(cbs_confidence - euk_confidence, 4)
  ) %>%
  relocate(rank, higher_rank, dataset,
           euk_cutoff, cbs_cutoff, diff_cutoff,
           euk_confidence, cbs_confidence, diff_confidence,
           euk_n_seq, cbs_n_seq, euk_n_group, cbs_n_group,
           euk_max_prop, cbs_max_prop)

cat(sprintf("Matched (rank, higher_rank, dataset) triples in both datasets: %d\n",
            nrow(joined)))
cat("\nMatches per (rank, higher_rank):\n")
joined %>% count(rank, higher_rank, sort = TRUE) %>% print(n = Inf)

dir.create("tests/CBSITS_eval/results", showWarnings = FALSE, recursive = TRUE)
write_tsv(joined, "tests/CBSITS_eval/results/full_comparison.tsv")

# ── Overall summary ─────────────────────────────────────────────────────────
cat("\n=== Overall summary (all matched triples) ===\n")
joined %>%
  summarise(
    n = n(),
    median_euk_conf = median(euk_confidence), median_cbs_conf = median(cbs_confidence),
    mean_diff_conf  = round(mean(diff_confidence), 4),
    pct_conf_improved = round(100 * mean(diff_confidence > 0), 1),
    mean_diff_cutoff = round(mean(diff_cutoff), 4),
    pct_cutoff_higher = round(100 * mean(diff_cutoff > 0), 1)
  ) %>%
  print()

cat("\n=== Summary by target rank ===\n")
joined %>%
  group_by(rank) %>%
  summarise(
    n = n(),
    median_euk_conf = round(median(euk_confidence), 3),
    median_cbs_conf = round(median(cbs_confidence), 3),
    mean_diff_conf  = round(mean(diff_confidence), 4),
    pct_conf_improved = round(100 * mean(diff_confidence > 0), 1),
    mean_diff_cutoff = round(mean(diff_cutoff), 4)
  ) %>%
  arrange(factor(rank, levels = c("species","genus","family","order","class","phylum"))) %>%
  print(n = Inf)

# ── Selection of ~30 groups spanning different levels ──────────────────────
# For each (rank, higher_rank) pair present in the intersection, keep the
# top few groups by CBS sequence count (best-sampled -> most trustworthy CBS
# estimate), so the final table is a reproducible, non-cherry-picked sample
# spanning species->genus all the way to class->phylum plus the two global
# (whole-database) rows.

level_order <- tibble(
  rank         = c("species","species","genus","genus","family","family","order","class","species","genus","family","order"),
  higher_rank  = c("genus","family","family","order","order","class","class","phylum","global","global","global","global"),
  level_rank   = 1:12
)

per_pair_n <- 4L

selected <- joined %>%
  inner_join(level_order, by = c("rank", "higher_rank")) %>%
  group_by(rank, higher_rank) %>%
  arrange(desc(cbs_n_seq), .by_group = TRUE) %>%
  slice_head(n = per_pair_n) %>%
  ungroup() %>%
  arrange(level_rank, desc(cbs_n_seq)) %>%
  select(-level_rank)

cat(sprintf("\n=== Selected comparison table: %d rows ===\n", nrow(selected)))
print(as.data.frame(selected), row.names = FALSE)

write_tsv(selected, "tests/CBSITS_eval/results/selected_30_comparison.tsv")
cat("\nWrote:\n  tests/CBSITS_eval/results/full_comparison.tsv\n  tests/CBSITS_eval/results/selected_30_comparison.tsv\n")
