#!/usr/bin/env Rscript
# test_script_params.R — Static checks that this repo's pipeline scripts carry
# the agreed similarity and sampling settings. Not synced to the sibling repos,
# whose scripts are their own.
#
# Usage: Rscript tests/test_script_params.R
# Note:  Must be run from the project root directory.

n_pass <- 0L
n_fail <- 0L
check <- function(desc, ok) {
  if (isTRUE(ok)) {
    n_pass <<- n_pass + 1L
    cat("PASS ", desc, "\n", sep = "")
  } else {
    n_fail <<- n_fail + 1L
    cat("FAIL ", desc, "\n", sep = "")
  }
}

# ── The shell scripts carry the same per-region values ────────────────────────
# 04, 06, 07a and 07c each declare their own MIN_OVERLAP, so a change to one
# must be made in all four. IDDEF must be 2 in each, and each must pass the
# overlap to its R script.

for (sc in c("scripts/04_remove_complexes.sh", "scripts/06_compute_sim.sh",
             "scripts/07a_predict_cutoffs.sh", "scripts/07c_predict_cutoffs_parallel.sh")) {
  lines <- readLines(sc)
  name  <- basename(sc)
  check(sprintf("%s: IDDEF=2", name), any(lines == "readonly IDDEF=2" | grepl("^readonly IDDEF=2( |$)", lines)))
  check(sprintf("%s: MIN_OVERLAP full_ITS=250, ITS1=50, ITS2=50", name),
        any(lines == 'MIN_OVERLAP["full_ITS"]=250') && any(lines == 'MIN_OVERLAP["ITS1"]=50') &&
        any(lines == 'MIN_OVERLAP["ITS2"]=50'))
  n_flag <- sum(grepl('--min_overlap +"\\$\\{MIN_OVERLAP\\[\\$region\\]\\}"', lines))
  check(sprintf("%s: passes --min_overlap per region (%d call site%s)", name, n_flag,
                if (n_flag == 1L) "" else "s"),
        n_flag == (if (grepl("07", name)) 2L else 1L))
}
for (sc in c("scripts/07a_predict_cutoffs.sh", "scripts/07c_predict_cutoffs_parallel.sh")) {
  check(sprintf("%s: does not pass --max_seq_no (predict.R's default applies)", basename(sc)),
        !any(grepl("--max_seq_no|MAX_SEQ_NO", readLines(sc))))
}
sub_lines <- readLines("scripts/05_prepare_subsets.sh")
check("05_prepare_subsets.sh: whole-species protocol (20000 sequences, max proportion 0.7, kingdom proportion 0.5, per-species cap 10, singleton cap 0.7)",
      any(sub_lines == "readonly MAX_SEQUENCES=20000") && any(sub_lines == "readonly MAX_PROPORTION=0.7") &&
      any(sub_lines == "readonly MAX_KINGDOM_PROPORTION=0.5") &&
      any(sub_lines == "readonly MAX_SEQS_PER_GROUP=10") && any(sub_lines == "readonly MAX_SINGLETON_PROPORTION=0.7"))
m_seq <- as.integer(sub("readonly MAX_SEQUENCES=", "", sub_lines[grepl("^readonly MAX_SEQUENCES=", sub_lines)]))
pred_src <- readLines("R/predict.R")
opt_at   <- grep('make_option("--max_seq_no"', pred_src, fixed = TRUE)
m_pred   <- as.integer(sub(".*default = ([0-9]+)L?.*", "\\1",
                           grep("default =", pred_src[opt_at + 0:2], value = TRUE)[1]))
check("predict.R's default --max_seq_no is at or above subset.R's --max_sequences (05), so it never triggers",
      length(m_seq) == 1L && length(m_pred) == 1L && m_pred >= m_seq)

# ── Result ────────────────────────────────────────────────────────────────────

cat(sprintf("\nRESULT: %d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
