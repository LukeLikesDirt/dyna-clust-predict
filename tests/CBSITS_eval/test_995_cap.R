#!/usr/bin/env Rscript
# test_995_cap.R — Test the user's other idea: hard-cap species-level cutoffs
# at 0.995. Uses the cached per-threshold F-measure traces already computed
# by the real production run (data/full_ITS/eukaryome.predicted, iddef=2,
# full-size max_seq_no=20000 data) -- no vsearch re-run needed, since the
# F-measure at every 0.001 threshold step from 0.90 to 1.00 is already stored
# for every dataset.
#
# For every species-level dataset currently at cutoff >= 0.995, this reports
# what the cutoff and F-measure would be if forced down to 0.995 instead of
# whatever threshold the sweep actually picked as optimal.
#
# This is a TEST ONLY. It does not modify R/predict.R or any pipeline script.
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

pred <- fromJSON("data/full_ITS/eukaryome.predicted", simplifyVector = FALSE)
sp <- pred$species

rows <- list()
for (name in names(sp)) {
  e <- sp[[name]]
  cutoff <- e[["cut-off"]]
  if (is.null(cutoff) || cutoff < 0.995) next

  fm <- e$fmeasures
  fkey <- "0.9950"
  f_at_995 <- if (!is.null(fm[[fkey]])) fm[[fkey]] else NA_real_

  rows[[length(rows) + 1]] <- data.table(
    higher_rank = e[["higher_rank"]],
    dataset = name,
    n_seq = e[["sequence number"]],
    n_grp = e[["group number"]],
    n_multiseq_grp = if (!is.null(e[["multiseq group number"]])) e[["multiseq group number"]] else NA_integer_,
    cut_actual = cutoff,
    f_actual = e[["confidence"]],
    cut_995 = 0.995,
    f_995 = f_at_995
  )
}
res <- rbindlist(rows, fill = TRUE)
res[, delta_f := f_995 - f_actual]
setorder(res, n_seq)

cat(sprintf("Species-level datasets with cutoff >= 0.995: %d\n\n", nrow(res)))
print(as.data.frame(res), row.names = FALSE)

cat("\n=== Overall impact of forcing cutoff -> 0.995 ===\n")
cat(sprintf("Rows affected (cutoff currently > 0.995, so a cap would change them): %d\n",
            sum(res$cut_actual > 0.995)))
cat(sprintf("Rows already exactly at cutoff == 0.995 (unaffected by cap): %d\n",
            sum(res$cut_actual == 0.995)))
cat(sprintf("Median F delta (995 - actual): %.4f\n", median(res$delta_f, na.rm = TRUE)))
cat(sprintf("Mean F delta   (995 - actual): %.4f\n", mean(res$delta_f, na.rm = TRUE)))
cat(sprintf("Rows where F IMPROVES at 0.995 (curve not monotonic near peak): %d / %d\n",
            sum(res$delta_f > 0, na.rm = TRUE), nrow(res)))
cat(sprintf("Rows where F drops > 0.05 at 0.995: %d / %d\n",
            sum(res$delta_f < -0.05, na.rm = TRUE), nrow(res)))
cat(sprintf("Rows where F drops > 0.10 at 0.995: %d / %d\n",
            sum(res$delta_f < -0.10, na.rm = TRUE), nrow(res)))

cat("\n=== By higher_rank ===\n")
print(res[, .(n = .N,
             median_delta_f = round(median(delta_f, na.rm = TRUE), 4),
             mean_delta_f = round(mean(delta_f, na.rm = TRUE), 4),
             n_big_loss_gt_0.05 = sum(delta_f < -0.05, na.rm = TRUE)),
        by = higher_rank][order(higher_rank)])

cat("\n=== Singleton-structure check: does multiseq group count explain who's ",
    "hurt worst by the cap? ===\n", sep = "")
res[, singleton_dominant := n_multiseq_grp <= 1 | is.na(n_multiseq_grp)]
print(res[, .(n = .N, median_delta_f = round(median(delta_f, na.rm = TRUE), 4)),
        by = singleton_dominant])

fwrite(res, "tests/CBSITS_eval/results/cap_995_test.tsv", sep = "\t")
cat("\nWrote tests/CBSITS_eval/results/cap_995_test.tsv\n")
