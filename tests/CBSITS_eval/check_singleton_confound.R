#!/usr/bin/env Rscript
# check_singleton_confound.R — Test whether the CBSITS confidence gain reported
# in tests/CBSITS_eval/report.html is genuine curation signal, or an artifact of
# CBSITS having a different singleton-species structure than EUKARYOME.
#
# Motivation: predict.R's Dice F-measure hands every singleton class a free
# perfect score at threshold 1.0, so datasets rich in single-sequence species
# score higher regardless of data quality. If CBSITS's matched datasets are more
# singleton-heavy than EUKARYOME's, part of the +0.17 median confidence gain is
# a scoring artifact rather than a curation effect.
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(data.table)
})

source("R/utils.R")

cmp <- read_tsv("tests/CBSITS_eval/results/full_comparison.tsv",
                show_col_types = FALSE)

euk_cls <- fread("data/ITS2/eukaryome_ITS2.classification", data.table = TRUE)
cbs_cls <- fread("data/CBSITS_ITS2/cbsits_ITS2.classification", data.table = TRUE)

parent_col <- c(genus = "genus", family = "family", order = "order",
                class = "class", phylum = "phylum", kingdom = "kingdom")

# Singleton fraction of the TARGET rank's classes within each parent taxon.
composition <- function(cls, target, hr) {
  col <- parent_col[[hr]]
  if (is.null(col) || !col %in% names(cls)) return(NULL)
  ok <- cls[is_identified(get(target)) & is_identified(get(col))]
  if (nrow(ok) == 0) return(NULL)
  counts <- ok[, .(k = .N), by = c(col, target)]
  setnames(counts, col, "parent")
  counts[, .(singleton_frac = mean(k == 1), seq_per_class = mean(k),
             n_class = .N), by = parent]
}

combos <- cmp %>% filter(higher_rank != "global") %>% distinct(rank, higher_rank)

enriched <- rbindlist(lapply(seq_len(nrow(combos)), function(i) {
  tgt <- combos$rank[i]; hr <- combos$higher_rank[i]
  e <- composition(euk_cls, tgt, hr); c <- composition(cbs_cls, tgt, hr)
  if (is.null(e) || is.null(c)) return(NULL)
  rows <- as.data.table(cmp)[rank == tgt & higher_rank == hr]
  rows <- merge(rows, e[, .(dataset = parent, euk_singleton = singleton_frac,
                            euk_seq_per_class = seq_per_class)], by = "dataset")
  merge(rows, c[, .(dataset = parent, cbs_singleton = singleton_frac,
                    cbs_seq_per_class = seq_per_class)], by = "dataset")
}), fill = TRUE)

cat("Matched datasets with composition data:", nrow(enriched), "\n\n")

cat("=== Singleton fraction: EUKARYOME vs CBSITS (same taxa) ===\n")
enriched %>%
  summarise(
    median_euk_singleton = round(median(euk_singleton), 3),
    median_cbs_singleton = round(median(cbs_singleton), 3),
    median_diff          = round(median(cbs_singleton - euk_singleton), 3),
    pct_cbs_more_singleton = round(100 * mean(cbs_singleton > euk_singleton), 1)
  ) %>% as.data.frame() %>% print(row.names = FALSE)

cat("\n=== Sequences per class ===\n")
enriched %>%
  summarise(
    median_euk_seq_per_class = round(median(euk_seq_per_class), 2),
    median_cbs_seq_per_class = round(median(cbs_seq_per_class), 2)
  ) %>% as.data.frame() %>% print(row.names = FALSE)

cat("\n=== Does the singleton shift explain the confidence gain? ===\n")
enriched[, singleton_shift := cbs_singleton - euk_singleton]
cat("Spearman(singleton_shift, diff_confidence) =",
    round(cor(enriched$singleton_shift, enriched$diff_confidence,
              method = "spearman"), 3), "\n")

fit <- lm(diff_confidence ~ singleton_shift, data = enriched)
cat("Linear fit: diff_confidence =",
    round(coef(fit)[1], 4), "+", round(coef(fit)[2], 4),
    "* singleton_shift   (R^2 =", round(summary(fit)$r.squared, 3), ")\n")
cat("  -> intercept is the confidence gain expected at ZERO singleton shift\n")

cat("\n=== Confidence gain in the subset where CBS is NOT more singleton-rich ===\n")
strict <- enriched[cbs_singleton <= euk_singleton]
cat("n =", nrow(strict), "\n")
strict %>%
  summarise(
    median_euk_conf = round(median(euk_confidence), 3),
    median_cbs_conf = round(median(cbs_confidence), 3),
    mean_diff_conf  = round(mean(diff_confidence), 4),
    pct_improved    = round(100 * mean(diff_confidence > 0), 1)
  ) %>% as.data.frame() %>% print(row.names = FALSE)

cat("\n=== Same, restricted to species-rank rows only ===\n")
sp <- strict[rank == "species"]
cat("n =", nrow(sp), "\n")
if (nrow(sp) > 0) {
  sp %>%
    summarise(
      median_euk_conf = round(median(euk_confidence), 3),
      median_cbs_conf = round(median(cbs_confidence), 3),
      mean_diff_conf  = round(mean(diff_confidence), 4),
      pct_improved    = round(100 * mean(diff_confidence > 0), 1)
    ) %>% as.data.frame() %>% print(row.names = FALSE)
}

write_tsv(enriched, "tests/CBSITS_eval/results/singleton_confound_check.tsv")
cat("\nWrote tests/CBSITS_eval/results/singleton_confound_check.tsv\n")
