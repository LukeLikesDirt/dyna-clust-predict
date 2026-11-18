#!/usr/bin/env Rscript
# test_singleton_exclusion.R — Prototype: what happens to the optimal cut-off if
# singleton classes are excluded from the F-measure sum but kept in the
# clustering graph?
#
# predict.R scores every class with a Dice term and weights it by class size:
#   dice(c) = max_r 2*count(r,c) / (class_size(c) + cluster_size(r))
# After dereplication every sequence is unique, so at threshold 1.0 each
# sequence is its own cluster and a class of size 1 scores 2*1/(1+1) = 1.0 --
# a free perfect score that no lower threshold can beat. Datasets dominated by
# single-sequence species therefore have their optimum pinned at 1.0 by
# construction, and their reported "confidence" is arithmetic, not signal.
#
# A singleton class says nothing about WHERE the within/between boundary lies,
# but it must stay in the graph because it can still wrongly absorb another
# species' sequences. So: keep singletons as distractors, drop them from the
# score.
#
# This script recomputes both scorings over the same similarity data and prints
# the two optima side by side.
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({
  library(data.table)
})

source("R/utils.R")

fasta_all <- "data/full_ITS/eukaryome_ITS.fasta"
class_all <- "data/full_ITS/eukaryome_ITS.classification"
tmp_dir   <- "tmp/singleton_test"
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

targets <- list(
  list(name = "Caladenia",    rank = "genus"),
  list(name = "Cypripedium",  rank = "genus"),
  list(name = "Pterostylis",  rank = "genus"),
  list(name = "Cladosporium", rank = "genus"),
  list(name = "Alternaria",   rank = "genus"),
  list(name = "Fusarium",     rank = "genus"),
  list(name = "Trichoderma",  rank = "genus")
)

start_t <- 0.8; end_t <- 1.0; step_t <- 0.001

cls <- fread(class_all)

# ── Sweep both scorings over one dataset ────────────────────────────────────
sweep_both <- function(sim, class_of, class_size, n) {
  # union-find over sequences, edges added as the threshold descends
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

  # scored: TRUE for classes counted in the F sum
  fm <- function(scored) {
    active <- which(is_root); tot <- 0; wt <- 0
    for (ci in seq_len(n_classes)) {
      if (class_size[ci] == 0L || !scored[ci]) next
      best <- 0
      for (r in active) {
        cnt <- root_cc[[r]][ci]
        if (cnt > 0L) {
          d <- 2 * cnt / (class_size[ci] + csize[r])
          if (d > best) best <- d
        }
      }
      tot <- tot + class_size[ci] * best
      wt  <- wt + class_size[ci]
    }
    if (wt == 0) return(NA_real_)
    tot / wt
  }

  scored_all   <- rep(TRUE, n_classes)
  scored_multi <- class_size >= 2L

  ord <- order(-sim$score)
  ii <- sim$i[ord]; jj <- sim$j[ord]; ss <- sim$score[ord]
  ptr <- 1L; ne <- length(ss)

  ts <- rev(seq(start_t, end_t, by = step_t))
  f_all <- numeric(length(ts)); f_multi <- numeric(length(ts))
  for (k in seq_along(ts)) {
    t <- ts[k]
    while (ptr <= ne && ss[ptr] >= t - 1e-9) { union(ii[ptr], jj[ptr]); ptr <- ptr + 1L }
    f_all[k]   <- fm(scored_all)
    f_multi[k] <- fm(scored_multi)
  }
  # back to ascending
  data.table(threshold = rev(ts), f_all = rev(f_all), f_multi = rev(f_multi))
}

pick <- function(ts, fs) {
  if (all(is.na(fs))) return(list(cut = NA, f = NA))
  best <- max(fs, na.rm = TRUE)
  tied <- which(fs == best)
  list(cut = ts[tied[ceiling(length(tied) / 2)]], f = best)
}

results <- list()

for (tg in targets) {
  sub <- cls[get(tg$rank) == tg$name & is_identified(species)]
  if (nrow(sub) < 30) { cat("skip", tg$name, "\n"); next }

  ids_file <- file.path(tmp_dir, paste0(tg$name, ".ids"))
  fa_file  <- file.path(tmp_dir, paste0(tg$name, ".fasta"))
  sim_file <- file.path(tmp_dir, paste0(tg$name, ".sim"))

  writeLines(sub$id, ids_file)
  system(sprintf("seqkit grep -f '%s' '%s' > '%s' 2>/dev/null",
                 ids_file, fasta_all, fa_file))
  system(sprintf(paste("vsearch --allpairs_global '%s' --acceptall --userout '%s'",
                       "--userfields query+target+id --threads 10"),
                 fa_file, sim_file), ignore.stdout = TRUE, ignore.stderr = TRUE)

  if (!file.exists(sim_file) || file.size(sim_file) == 0) { cat("no sim for", tg$name, "\n"); next }

  raw <- fread(sim_file, header = FALSE, col.names = c("q", "t", "pid"))
  raw[, score := pid / 100]

  ids <- unique(c(raw$q, raw$t, sub$id))
  idx <- setNames(seq_along(ids), ids)
  sp  <- setNames(sub$species, sub$id)[ids]
  classes <- unique(sp[!is.na(sp)])
  class_of <- match(sp, classes); class_of[is.na(class_of)] <- 0L
  class_size <- tabulate(class_of[class_of > 0], nbins = length(classes))

  sim <- data.table(i = idx[raw$q], j = idx[raw$t], score = raw$score)

  tr <- sweep_both(sim, class_of, class_size, length(ids))
  a <- pick(tr$threshold, tr$f_all)
  m <- pick(tr$threshold, tr$f_multi)

  results[[tg$name]] <- data.table(
    dataset        = tg$name,
    n_seq          = length(ids),
    n_species      = length(classes),
    singleton_frac = round(mean(class_size == 1), 3),
    cut_current    = a$cut, f_current = round(a$f, 4),
    cut_nosingle   = m$cut, f_nosingle = round(m$f, 4)
  )
  cat(sprintf("%-14s current: cut=%.3f F=%.4f   |   singletons excluded: cut=%.3f F=%.4f\n",
              tg$name, a$cut, a$f, m$cut, m$f))
}

out <- rbindlist(results)
cat("\n=== Summary ===\n")
print(as.data.frame(out), row.names = FALSE)

dir.create("tests/CBSITS_eval/results", showWarnings = FALSE, recursive = TRUE)
fwrite(out, "tests/CBSITS_eval/results/singleton_exclusion_test.tsv", sep = "\t")
cat("\nWrote tests/CBSITS_eval/results/singleton_exclusion_test.tsv\n")
