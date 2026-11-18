#!/usr/bin/env Rscript
# test_iddef_pinned_species.R — For every species-level cutoff=1 row in the
# CURRENT (post-completeness-filter, post-complex-removal) full_ITS cutoffs
# table, compare --iddef 1 vs --iddef 2 (the pipeline's current default,
# already baked into data/full_ITS/eukaryome.cutoffs.json.txt).
#
# Revisits the iddef question from earlier in this project under materially
# different conditions: complex removal now strips genuinely indistinguishable
# species pairs out of the dataset before prediction ever runs, and the
# completeness filter + min_multiseq_groups guard remove the worst truncation/
# singleton artifacts. The earlier iddef=1 vs iddef=2 comparison (which found
# iddef=1 costs 0.05-0.15 F-measure) was run on data with none of those fixes
# applied -- worth checking whether that cost still holds now that a lot of
# the noise it may have been (correctly) reacting to has already been removed
# by other means.
#
# Input: data/full_ITS/eukaryome.cutoffs.json.txt (for the list of pinned
# groups) and data/full_ITS/eukaryome_ITS_nocomplex.{fasta,classification}
# (the actual post-fix sequence data to re-test on).
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({
  library(data.table)
})

source("R/utils.R")

fasta_all <- "data/full_ITS/eukaryome_ITS_nocomplex.fasta"
class_all <- "data/full_ITS/eukaryome_ITS_nocomplex.classification"
cutoffs_in <- "data/full_ITS/eukaryome.cutoffs.json.txt"
tmp_dir   <- "tmp/iddef_pinned_test"
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

start_t <- 0.90; end_t <- 1.0; step_t <- 0.001
# Capped at 5000 (not predict.R's 20000) for tractability running locally on
# 12 cores: allpairs_global is O(n^2), and a full run of all 72 groups at
# their production sizes (up to 20000, several already capped there) was
# timed at ~6.3hr for both iddef settings combined. At 5000 it's ~70min and
# still a representative random sample of the same pool predict.R itself drew
# its (also-capped) sample from -- only the largest 7 of 72 groups (the two
# phyla, kingdom Fungi, and 4 of the 6 classes/orders at >5000) are affected.
MAX_SEQ_NO <- 5000L

parent_col <- c(genus = "genus", family = "family", order = "order",
                class = "class", phylum = "phylum", kingdom = "kingdom")

# ── Load the list of pinned (higher_rank, dataset) groups ──────────────────

cutoffs <- fread(cutoffs_in, sep = "\t", header = TRUE)
pinned <- cutoffs[rank == "species" & `cut-off` == 1,
                  .(higher_rank, dataset, orig_n_seq = `sequence number`,
                    orig_confidence = confidence)]
setorder(pinned, orig_n_seq)
cat(sprintf("Found %d species-level cutoff=1 groups to re-test.\n", nrow(pinned)))
print(as.data.frame(pinned))

cls <- fread(class_all, sep = "\t", header = TRUE, colClasses = "character")

# ── Union-find F-measure sweep (same logic used throughout this project's
#    iddef investigation) ────────────────────────────────────────────────────

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
        if (cnt > 0L) { d <- 2 * cnt / (class_size[ci] + csize[r]); if (d > best) best <- d }
      }
      tot <- tot + class_size[ci] * best
    }
    tot / sum(class_size)
  }
  if (nrow(sim) == 0) return(list(cut = NA, f = NA))
  ord <- order(-sim$score); ii <- sim$i[ord]; jj <- sim$j[ord]; ss <- sim$score[ord]
  ptr <- 1L; ne <- length(ss)
  ts <- rev(seq(start_t, end_t, by = step_t)); fs <- numeric(length(ts))
  for (k in seq_along(ts)) {
    t <- ts[k]
    while (ptr <= ne && ss[ptr] >= t - 1e-9) { union(ii[ptr], jj[ptr]); ptr <- ptr + 1L }
    fs[k] <- fm()
  }
  ts <- rev(ts); fs <- rev(fs)
  best <- max(fs); tied <- which(fs == best)
  list(cut = ts[tied[ceiling(length(tied) / 2)]], f = round(best, 4))
}

run_one_iddef <- function(fa_file, sim_file, ids, sp_map, iddef, n_threads) {
  system(sprintf(paste("vsearch --allpairs_global '%s' --acceptall --userout '%s'",
                       "--userfields query+target+id --iddef %d --threads %d"),
                 fa_file, sim_file, iddef, n_threads),
         ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (!file.exists(sim_file) || file.size(sim_file) == 0) return(list(cut = NA, f = NA, n_seq = length(ids)))
  raw <- fread(sim_file, header = FALSE, col.names = c("q", "t", "pid"),
              colClasses = c("character", "character", "numeric"))
  raw[, score := pid / 100]
  idx <- setNames(seq_along(ids), ids)
  sp <- sp_map[ids]
  classes <- unique(sp[!is.na(sp)])
  class_of <- match(sp, classes); class_of[is.na(class_of)] <- 0L
  class_size <- tabulate(class_of[class_of > 0], nbins = length(classes))
  res <- sweep_f(data.table(i = idx[raw$q], j = idx[raw$t], score = raw$score),
                 class_of, class_size, length(ids))
  c(res, list(n_seq = length(ids), n_spp = length(classes)))
}

# ── Process each pinned group ────────────────────────────────────────────────

out <- list()
for (i in seq_len(nrow(pinned))) {
  hr <- pinned$higher_rank[i]; ds <- pinned$dataset[i]
  col <- parent_col[[hr]]
  sub <- cls[get(col) == ds & is_identified(species)]
  if (nrow(sub) < 10) { cat("skip (too few):", hr, ds, "\n"); next }

  ids <- sub$id
  if (length(ids) > MAX_SEQ_NO) {
    set.seed(42)
    ids <- sample(ids, MAX_SEQ_NO)
  }
  sp_map <- setNames(sub$species, sub$id)

  safe_name <- gsub("[^A-Za-z0-9_]", "_", paste(hr, ds, sep = "_"))
  fa_file <- file.path(tmp_dir, paste0(safe_name, ".fasta"))
  n_threads <- if (length(ids) > 5000) 10L else if (length(ids) > 500) 6L else 2L

  writeLines(ids, file.path(tmp_dir, paste0(safe_name, ".ids")))
  system(sprintf("seqkit grep -f '%s' '%s' > '%s' 2>/dev/null",
                 file.path(tmp_dir, paste0(safe_name, ".ids")), fasta_all, fa_file))

  t0 <- Sys.time()
  r2 <- run_one_iddef(fa_file, file.path(tmp_dir, paste0(safe_name, "_id2.sim")),
                      ids, sp_map, 2L, n_threads)
  r1 <- run_one_iddef(fa_file, file.path(tmp_dir, paste0(safe_name, "_id1.sim")),
                      ids, sp_map, 1L, n_threads)
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)

  unlink(fa_file); unlink(file.path(tmp_dir, paste0(safe_name, ".ids")))
  unlink(file.path(tmp_dir, paste0(safe_name, "_id2.sim")))
  unlink(file.path(tmp_dir, paste0(safe_name, "_id1.sim")))

  cat(sprintf("[%2d/%2d] %-8s %-22s n=%-6d (%.1fs) | iddef2: cut=%s F=%s | iddef1: cut=%s F=%s\n",
              i, nrow(pinned), hr, ds, r2$n_seq, elapsed,
              ifelse(is.na(r2$cut), "NA", sprintf("%.3f", r2$cut)), ifelse(is.na(r2$f), "NA", sprintf("%.4f", r2$f)),
              ifelse(is.na(r1$cut), "NA", sprintf("%.3f", r1$cut)), ifelse(is.na(r1$f), "NA", sprintf("%.4f", r1$f))))

  out[[length(out) + 1]] <- data.table(
    higher_rank = hr, dataset = ds, n_seq = r2$n_seq, n_spp = r2$n_spp,
    cut_iddef2 = r2$cut, f_iddef2 = r2$f,
    cut_iddef1 = r1$cut, f_iddef1 = r1$f,
    elapsed_sec = elapsed
  )

  # Write incrementally so partial progress survives if the run is interrupted
  fwrite(rbindlist(out), "tests/CBSITS_eval/results/iddef_pinned_species.tsv", sep = "\t")
}

res <- rbindlist(out)
cat("\n=== SUMMARY ===\n")
print(as.data.frame(res), row.names = FALSE)

cat("\n=== How many un-pin from 1.0 under iddef=1? ===\n")
cat(sprintf("Still at cut=1.000 under iddef=1: %d / %d\n",
            sum(res$cut_iddef1 >= 0.9995, na.rm = TRUE), nrow(res)))

cat("\n=== F-measure delta (iddef1 - iddef2) ===\n")
res[, delta_f := f_iddef1 - f_iddef2]
print(res[, .(median_delta_f = round(median(delta_f, na.rm = TRUE), 4),
             mean_delta_f = round(mean(delta_f, na.rm = TRUE), 4),
             pct_iddef1_better = round(100 * mean(delta_f > 0, na.rm = TRUE), 1)),
        by = higher_rank])

cat("\nWrote tests/CBSITS_eval/results/iddef_pinned_species.tsv\n")
