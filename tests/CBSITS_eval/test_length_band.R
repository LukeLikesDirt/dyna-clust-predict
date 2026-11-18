#!/usr/bin/env Rscript
# test_length_band.R — Does restricting full-ITS records to a narrow length band
# rescue the species-level cut-off from being pinned at 1.0?
#
# vsearch --allpairs_global does global (end-to-end) alignment and, unlike
# usearch_global, accepts no --strand or coverage options. EUKARYOME's full-ITS
# records run 250-1500 bp, so a partial record aligned against a complete one of
# the SAME species is heavily penalised for the overhang. Measured on Fusarium:
# conspecific pairs whose length ratio is < 0.6 have a median identity of 44%,
# against 89% for pairs of comparable length. Conspecifics that cannot exceed
# the sweep's lower bound never cluster, so "every sequence its own species"
# (cut-off 1.0) wins by default.
#
# This restricts each genus to the core of its own length distribution and
# re-runs the sweep, reporting cut-off and F-measure before and after.
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({
  library(data.table)
})

source("R/utils.R")

fasta_all <- "data/full_ITS/eukaryome_ITS.fasta"
class_all <- "data/full_ITS/eukaryome_ITS.classification"
tmp_dir   <- "tmp/lenband_test"
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

genera  <- c("Fusarium", "Alternaria", "Trichoderma", "Cladosporium")
start_t <- 0.8; end_t <- 1.0; step_t <- 0.001
# keep the central mass of each genus's own length distribution
lo_q <- 0.25; hi_q <- 0.75

cls <- fread(class_all)

sweep_f <- function(sim, class_of, class_size, n) {
  parent <- seq_len(n); rnk <- integer(n); csize <- rep(1L, n)
  n_classes <- length(class_size)
  root_cc <- lapply(seq_len(n), function(i) {
    v <- integer(n_classes); if (class_of[i] > 0L) v[class_of[i]] <- 1L; v
  })
  is_root <- rep(TRUE, n)
  find <- function(x) { while (parent[x] != x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }; x }
  union <- function(a, b) {
    ra <- find(a); rb <- find(b); if (ra == rb) return(invisible())
    if (rnk[ra] < rnk[rb]) { t <- ra; ra <- rb; rb <- t }
    parent[rb] <<- ra
    if (rnk[ra] == rnk[rb]) rnk[ra] <<- rnk[ra] + 1L
    csize[ra] <<- csize[ra] + csize[rb]
    root_cc[[ra]] <<- root_cc[[ra]] + root_cc[[rb]]
    is_root[rb] <<- FALSE
  }
  fm <- function() {
    active <- which(is_root); tot <- 0
    for (ci in seq_len(n_classes)) {
      if (class_size[ci] == 0L) next
      best <- 0
      for (r in active) {
        cnt <- root_cc[[r]][ci]
        if (cnt > 0L) {
          d <- 2 * cnt / (class_size[ci] + csize[r]); if (d > best) best <- d
        }
      }
      tot <- tot + class_size[ci] * best
    }
    tot / sum(class_size)
  }
  ord <- order(-sim$score)
  ii <- sim$i[ord]; jj <- sim$j[ord]; ss <- sim$score[ord]
  ptr <- 1L; ne <- length(ss)
  ts <- rev(seq(start_t, end_t, by = step_t)); fs <- numeric(length(ts))
  for (k in seq_along(ts)) {
    t <- ts[k]
    while (ptr <= ne && ss[ptr] >= t - 1e-9) { union(ii[ptr], jj[ptr]); ptr <- ptr + 1L }
    fs[k] <- fm()
  }
  ts <- rev(ts); fs <- rev(fs)
  best <- max(fs); tied <- which(fs == best)
  list(cut = ts[tied[ceiling(length(tied) / 2)]], f = best)
}

run_one <- function(ids, tag, species_map) {
  ids_file <- file.path(tmp_dir, paste0(tag, ".ids"))
  fa_file  <- file.path(tmp_dir, paste0(tag, ".fasta"))
  sim_file <- file.path(tmp_dir, paste0(tag, ".sim"))
  writeLines(ids, ids_file)
  system(sprintf("seqkit grep -f '%s' '%s' > '%s' 2>/dev/null", ids_file, fasta_all, fa_file))
  system(sprintf(paste("vsearch --allpairs_global '%s' --acceptall --userout '%s'",
                       "--userfields query+target+id --threads 10"), fa_file, sim_file),
         ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (!file.exists(sim_file) || file.size(sim_file) == 0) return(NULL)
  raw <- fread(sim_file, header = FALSE, col.names = c("q", "t", "pid"))
  raw[, score := pid / 100]
  all_ids <- unique(c(raw$q, raw$t, ids))
  idx <- setNames(seq_along(all_ids), all_ids)
  sp <- species_map[all_ids]
  classes <- unique(sp[!is.na(sp)])
  class_of <- match(sp, classes); class_of[is.na(class_of)] <- 0L
  class_size <- tabulate(class_of[class_of > 0], nbins = length(classes))
  res <- sweep_f(data.table(i = idx[raw$q], j = idx[raw$t], score = raw$score),
                 class_of, class_size, length(all_ids))
  c(res, list(n_seq = length(all_ids), n_spp = length(classes)))
}

out <- list()
for (g in genera) {
  sub <- cls[genus == g & is_identified(species)]
  if (nrow(sub) < 30) next
  sp_map <- setNames(sub$species, sub$id)

  writeLines(sub$id, file.path(tmp_dir, "all.ids"))
  system(sprintf("seqkit grep -f '%s' '%s' > '%s' 2>/dev/null",
                 file.path(tmp_dir, "all.ids"), fasta_all, file.path(tmp_dir, "all.fasta")))
  lens <- fread(cmd = sprintf("seqkit fx2tab -nil '%s' 2>/dev/null", file.path(tmp_dir, "all.fasta")),
                header = FALSE)
  setnames(lens, c("id", "len")); lens[, id := sub(" .*", "", id)]

  qs <- quantile(lens$len, c(lo_q, hi_q))
  keep <- lens[len >= qs[1] & len <= qs[2], id]

  full <- run_one(sub$id, paste0(g, "_full"), sp_map)
  band <- run_one(keep,   paste0(g, "_band"), sp_map)
  if (is.null(full) || is.null(band)) next

  cat(sprintf("%-13s all lengths: n=%-5d cut=%.3f F=%.4f  |  band %d-%dbp: n=%-5d cut=%.3f F=%.4f\n",
              g, full$n_seq, full$cut, full$f, qs[1], qs[2], band$n_seq, band$cut, band$f))
  out[[g]] <- data.table(genus = g,
                         n_all = full$n_seq, cut_all = full$cut, f_all = round(full$f, 4),
                         len_lo = qs[1], len_hi = qs[2],
                         n_band = band$n_seq, cut_band = band$cut, f_band = round(band$f, 4))
}

res <- rbindlist(out)
cat("\n=== Summary ===\n"); print(as.data.frame(res), row.names = FALSE)
dir.create("tests/CBSITS_eval/results", showWarnings = FALSE, recursive = TRUE)
fwrite(res, "tests/CBSITS_eval/results/length_band_test.tsv", sep = "\t")
cat("\nWrote tests/CBSITS_eval/results/length_band_test.tsv\n")
