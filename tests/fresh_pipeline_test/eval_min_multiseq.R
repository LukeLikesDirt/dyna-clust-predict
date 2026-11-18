#!/usr/bin/env Rscript
# eval_min_multiseq.R — What does predict.R's --min_multiseq_groups (07a/07c
# MIN_MULTISEQ_GROUPS) add on top of 05's MAX_SINGLETON_PROPORTION?
#
# Reads the eukaryome.predicted files of a fresh_test.R stage 07 run (they hold
# every dataset, including those the gate withheld) and the heavy-stage table,
# then compares gate on (production, 10) with gate off (0):
#   1. raw cutoffs of the withheld datasets vs the reported ones;
#   2. the 08 step (R/consolidate_cutoffs.R) run on both tables. Gate on, the
#      withheld datasets are listed with confidence NA, so 08 enumerates them
#      but must fill them from an ancestor -- what production does for a taxon
#      with no reported value of its own.
# predict.R's other report filters (>= 10 groups, >= 30 sequences) apply to
# both arms. data/ is only read.
#
# Usage (from the project root):
#   Rscript tests/fresh_pipeline_test/eval_min_multiseq.R --work <fresh_test work dir>

suppressPackageStartupMessages({
  library(optparse); library(data.table); library(jsonlite)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--work", default = "tmp/fresh_pipeline_test"),
  make_option("--region", default = "full_ITS"),
  make_option("--results", default = "tests/fresh_pipeline_test/results"),
  make_option("--gate", type = "integer", default = 10L)
)))
`%||%` <- function(a, b) if (is.null(a)) b else a
region <- opt$region; gate <- opt$gate
class_file <- c(full_ITS = "data/full_ITS/eukaryome_ITS_nocomplex.classification",
                ITS1 = "data/ITS1/eukaryome_ITS1_nocomplex.classification",
                ITS2 = "data/ITS2/eukaryome_ITS2_nocomplex.classification")[[region]]
out_dir <- file.path(opt$work, "min_multiseq_eval"); dir.create(out_dir, showWarnings = FALSE)

# ── 1. Every dataset predict.R computed ──────────────────────────────────────
pred_root <- file.path(opt$work, "pred", region)
rows <- list()
for (combo_dir in list.dirs(pred_root, recursive = FALSE)) {
  pj <- file.path(combo_dir, "eukaryome.predicted")
  if (!file.exists(pj)) next
  combo <- sub("_", ":", basename(combo_dir))
  pred <- fromJSON(pj, simplifyVector = FALSE)
  for (rk in names(pred)) for (ds in names(pred[[rk]])) {
    e <- pred[[rk]][[ds]]
    rows[[length(rows) + 1]] <- data.table(
      combo = combo, rank = rk, higher_rank = e[["higher_rank"]] %||% "global", dataset = ds,
      cutoff = as.numeric(e[["cut-off"]] %||% NA), confidence = as.numeric(e[["confidence"]] %||% NA),
      seq_n = as.integer(e[["sequence number"]] %||% NA), grp_n = as.integer(e[["group number"]] %||% NA),
      multiseq_n = as.integer(e[["multiseq group number"]] %||% NA),
      max_prop = as.numeric(e[["max proportion"]] %||% NA))
  }
}
ds <- rbindlist(rows)
if (nrow(ds) == 0) stop("No eukaryome.predicted files under ", pred_root)

# Kingdom and global species datasets too large for predict.R (heavy stage)
heavy_file <- file.path(opt$results, sprintf("heavy_%s.tsv", region))
if (file.exists(heavy_file)) {
  h <- fread(heavy_file)
  ds <- rbind(ds, h[, .(combo = ifelse(dataset == "global", "species:global", "species:kng"),
                        rank = "species",
                        higher_rank = ifelse(dataset == "global", "global", "kingdom"),
                        dataset = ifelse(dataset == "global", "All", sub("@kng$", "", dataset)),
                        cutoff, confidence = f_at_cutoff, seq_n = n, grp_n = n_spp,
                        multiseq_n = n_multiseq, max_prop)])
}

ds <- ds[grp_n >= 10 & seq_n >= 30]                 # predict.R's other filters
ds[, gated := multiseq_n < gate]
ds[, singleton_pct := round(100 * (1 - multiseq_n / grp_n), 1)]
ds[, f_is_1 := confidence >= 0.99995]
fwrite(ds, file.path(opt$results, sprintf("min_multiseq_datasets_%s.tsv", region)), sep = "\t")

q <- function(x, p) { x <- x[!is.na(x)]; if (length(x)) as.numeric(quantile(x, p, names = FALSE)) else NA_real_ }
raw_summary <- ds[, .(datasets = .N, withheld = sum(gated),
  withheld_multiseq_0 = sum(gated & multiseq_n == 0),
  withheld_multiseq_1_9 = sum(gated & multiseq_n > 0),
  reported_median_cutoff = q(cutoff[!gated], 0.5),
  withheld_median_cutoff = q(cutoff[gated], 0.5),
  reported_median_F = q(confidence[!gated], 0.5),
  withheld_median_F = q(confidence[gated], 0.5),
  reported_pct_F1 = round(100 * mean(f_is_1[!gated]), 1),
  withheld_pct_F1 = round(100 * mean(f_is_1[gated]), 1),
  withheld_median_singleton_pct = q(singleton_pct[gated], 0.5)), by = combo]
fwrite(raw_summary, file.path(opt$results, sprintf("min_multiseq_raw_summary_%s.tsv", region)), sep = "\t")

bins <- ds[, .(datasets = .N, median_cutoff = q(cutoff, 0.5), median_F = q(confidence, 0.5),
               pct_F1 = round(100 * mean(f_is_1), 1), median_singleton_pct = q(singleton_pct, 0.5)),
           by = .(multiseq = cut(multiseq_n, c(-1, 0, 4, 9, 19, Inf),
                                 labels = c("0", "1-4", "5-9", "10-19", "20+")))][order(multiseq)]
fwrite(bins, file.path(opt$results, sprintf("min_multiseq_by_bin_%s.tsv", region)), sep = "\t")

# ── 2. Through the 08 step, both arms ────────────────────────────────────────
write_raw <- function(d, path) {
  out <- d[, .(rank, higher_rank, dataset, `cut-off` = cutoff, confidence,
               `sequence number` = seq_n, `group number` = grp_n,
               `multiseq group number` = multiseq_n, `max proportion` = max_prop)]
  fwrite(out, path, sep = "\t", na = "NA")
}
arm_on  <- copy(ds)[gated == TRUE, confidence := NA_real_]
arm_off <- copy(ds)
consolidate <- function(d, tag) {
  raw <- file.path(out_dir, sprintf("raw_%s.txt", tag)); res <- file.path(out_dir, sprintf("consolidated_%s.txt", tag))
  write_raw(d, raw)
  st <- system2("Rscript", c("R/consolidate_cutoffs.R", "--cutoffs_in", raw,
                             "--classification_in", class_file, "--output", res),
                stdout = file.path(out_dir, sprintf("consolidate_%s.log", tag)), stderr = file.path(out_dir, sprintf("consolidate_%s.log", tag)))
  if (st != 0) stop("consolidate_cutoffs.R failed (", tag, "); see ", out_dir)
  fread(res, sep = "\t", na.strings = "")
}
c_on  <- consolidate(arm_on,  "gate_on")
c_off <- consolidate(arm_off, "gate_off")
key <- c("rank", "higher_rank", "dataset")
cmp <- merge(c_on[,  .(rank, higher_rank, dataset, cutoff_on = `cut-off`, conf_on = confidence, source_on = source, clamped_on = clamped)],
             c_off[, .(rank, higher_rank, dataset, cutoff_off = `cut-off`, conf_off = confidence, source_off = source, clamped_off = clamped)],
             by = key, all = TRUE)
cmp <- merge(cmp, ds[, .(rank, higher_rank, dataset, own_withheld = gated, own_multiseq = multiseq_n)], by = key, all.x = TRUE)
cmp[is.na(own_withheld), own_withheld := FALSE]
cmp[, diff := round(cutoff_off - cutoff_on, 4)]
cmp[, cell := fifelse(own_withheld, "own value withheld by gate", "own value reported (or none computed)")]
fwrite(cmp, file.path(opt$results, sprintf("min_multiseq_consolidated_%s.tsv", region)), sep = "\t")

final_summary <- cmp[, .(cells = .N,
  changed = sum(abs(diff) > 1e-9, na.rm = TRUE),
  filled_only_gate_off = sum(is.na(cutoff_on) & !is.na(cutoff_off)),
  median_diff_changed = q(diff[abs(diff) > 1e-9], 0.5),
  q10_diff_changed = q(diff[abs(diff) > 1e-9], 0.1),
  q90_diff_changed = q(diff[abs(diff) > 1e-9], 0.9),
  max_abs_diff = suppressWarnings(max(abs(diff), na.rm = TRUE)),
  gate_off_self_wins = sum(source_off == "self", na.rm = TRUE),
  gate_off_self_wins_F1 = sum(source_off == "self" & conf_off >= 0.99995, na.rm = TRUE)),
  by = .(rank, higher_rank, cell)][order(rank, higher_rank, cell)]
fwrite(final_summary, file.path(opt$results, sprintf("min_multiseq_final_summary_%s.tsv", region)), sep = "\t")

options(width = 200)
cat("\n== Raw: datasets withheld by --min_multiseq_groups", gate, "(other filters applied to both) ==\n")
print(raw_summary)
cat("\n== Raw: by number of multi-sequence groups ==\n"); print(bins)
cat("\n== After 08 (consolidate_cutoffs.R): gate off minus gate on ==\n"); print(final_summary)
cat("\nWritten to", opt$results, "and", out_dir, "\n")
