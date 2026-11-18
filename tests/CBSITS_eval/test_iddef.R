#!/usr/bin/env Rscript
# test_iddef.R — Phase 1a+1b combined: does vsearch's --iddef (identity
# definition) explain pinned/inflated cut-offs, and does restricting to
# "complete ITS span" records compound with it?
#
# Background: compute_sim() in predict.R/compute_sim.R runs vsearch with no
# --iddef, so it defaults to 2 (edit distance EXCLUDING terminal gaps). A
# partial record that is a perfect substring of a complete one then scores
# 100% identity. Measured on Fusarium: --iddef 2 manufactures 995 "identical"
# cross-species pairs at 100%; --iddef 1 (terminal gaps COUNTED) finds zero.
#
# Pairwise vsearch scores for (i,j) do not depend on what else is aligned in
# the same run, so ONE sim matrix per (region, genus, iddef) suffices -- the
# "complete span only" comparison is done by subsetting that same matrix to
# IDs in the complete-span set, not by re-running vsearch.
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({
  library(data.table)
})

source("R/utils.R")

tmp_dir <- "tmp/iddef_test"
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

targets <- list(
  list(genus = "Fusarium",     region = "full_ITS"),
  list(genus = "Alternaria",   region = "full_ITS"),
  list(genus = "Trichoderma",  region = "full_ITS"),
  list(genus = "Cladosporium", region = "full_ITS"),
  list(genus = "Penicillium",  region = "full_ITS"),
  list(genus = "Aspergillus",  region = "full_ITS"),
  list(genus = "Caladenia",    region = "full_ITS"),
  list(genus = "Cypripedium",  region = "full_ITS"),
  list(genus = "Fusarium",     region = "ITS2"),
  list(genus = "Alternaria",   region = "ITS2"),
  list(genus = "Trichoderma",  region = "ITS2"),
  list(genus = "Cladosporium", region = "ITS2"),
  list(genus = "Penicillium",  region = "ITS2"),
  list(genus = "Aspergillus",  region = "ITS2"),
  list(genus = "Caladenia",    region = "ITS2"),
  list(genus = "Cypripedium",  region = "ITS2")
)

iddefs <- c(1L, 2L, 4L)
start_t <- 0.8; end_t <- 1.0; step_t <- 0.001

region_paths <- list(
  full_ITS = list(fasta = "data/full_ITS/eukaryome_ITS.fasta",
                  cls   = "data/full_ITS/eukaryome_ITS.classification"),
  ITS2     = list(fasta = "data/ITS2/eukaryome_ITS2.fasta",
                  cls   = "data/ITS2/eukaryome_ITS2.classification")
)

# ── Complete-span ID set: present in BOTH ITS1 and ITS2 outputs ────────────
cat("Building complete-span ID set (ITS1 ids intersect ITS2 ids)...\n")
get_ids <- function(f) sub("^>", "", grep("^>", readLines(f), value = TRUE))
its1_ids <- get_ids("data/ITS1/eukaryome_ITS1.fasta")
its2_ids <- get_ids("data/ITS2/eukaryome_ITS2.fasta")
complete_ids <- intersect(its1_ids, its2_ids)
cat("Complete-span IDs:", length(complete_ids), "\n")
complete_set <- new.env()
for (id in complete_ids) assign(id, TRUE, envir = complete_set)
is_complete <- function(ids) vapply(ids, function(x) exists(x, envir = complete_set, inherits = FALSE), logical(1))

# ── F-measure sweep over a similarity table, restricted to a subset of ids ──
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
  list(cut = ts[tied[ceiling(length(tied) / 2)]], f = round(best, 4))
}

run_sweep_on_subset <- function(raw_sim, sp_map, subset_ids) {
  keep_ids <- intersect(names(sp_map)[sp_map %in% sp_map], subset_ids)
  sim <- raw_sim[q %in% subset_ids & t %in% subset_ids]
  all_ids <- intersect(subset_ids, names(sp_map))
  if (length(all_ids) < 3 || nrow(sim) == 0) return(list(cut = NA, f = NA, n_seq = length(all_ids)))
  idx <- setNames(seq_along(all_ids), all_ids)
  sp <- sp_map[all_ids]
  classes <- unique(sp[!is.na(sp)])
  class_of <- match(sp, classes); class_of[is.na(class_of)] <- 0L
  class_size <- tabulate(class_of[class_of > 0], nbins = length(classes))
  sim2 <- sim[q %in% all_ids & t %in% all_ids]
  res <- sweep_f(data.table(i = idx[sim2$q], j = idx[sim2$t], score = sim2$score),
                 class_of, class_size, length(all_ids))
  c(res, list(n_seq = length(all_ids), n_spp = length(classes)))
}

out <- list()

for (tg in targets) {
  paths <- region_paths[[tg$region]]
  cls <- fread(paths$cls)
  sub <- cls[genus == tg$genus & is_identified(species)]
  if (nrow(sub) < 30) { cat("skip (too few seqs):", tg$genus, tg$region, "\n"); next }

  sp_map <- setNames(sub$species, sub$id)
  ids_file <- file.path(tmp_dir, paste0(tg$genus, "_", tg$region, ".ids"))
  fa_file  <- file.path(tmp_dir, paste0(tg$genus, "_", tg$region, ".fasta"))
  writeLines(sub$id, ids_file)
  system(sprintf("seqkit grep -f '%s' '%s' > '%s' 2>/dev/null", ids_file, paths$fasta, fa_file))

  comp_ids <- sub$id[sub$id %in% complete_ids]

  for (idd in iddefs) {
    sim_file <- file.path(tmp_dir, paste0(tg$genus, "_", tg$region, "_iddef", idd, ".sim"))
    system(sprintf(paste("vsearch --allpairs_global '%s' --acceptall --userout '%s'",
                         "--userfields query+target+id --iddef %d --threads 10"),
                   fa_file, sim_file, idd), ignore.stdout = TRUE, ignore.stderr = TRUE)
    if (!file.exists(sim_file) || file.size(sim_file) == 0) {
      cat("no sim:", tg$genus, tg$region, "iddef", idd, "\n"); next
    }
    raw <- fread(sim_file, header = FALSE, col.names = c("q", "t", "pid"))
    raw[, score := pid / 100]

    r_all  <- run_sweep_on_subset(raw, sp_map, sub$id)
    r_comp <- run_sweep_on_subset(raw, sp_map, comp_ids)

    cat(sprintf("%-13s %-9s iddef=%d | ALL n=%-5d cut=%s F=%s  |  COMPLETE n=%-5d cut=%s F=%s\n",
                tg$genus, tg$region, idd,
                r_all$n_seq,  ifelse(is.na(r_all$cut),  "NA", sprintf("%.3f", r_all$cut)),  ifelse(is.na(r_all$f),  "NA", sprintf("%.4f", r_all$f)),
                r_comp$n_seq, ifelse(is.na(r_comp$cut), "NA", sprintf("%.3f", r_comp$cut)), ifelse(is.na(r_comp$f), "NA", sprintf("%.4f", r_comp$f))))

    out[[length(out) + 1]] <- data.table(
      genus = tg$genus, region = tg$region, iddef = idd,
      n_all = r_all$n_seq, cut_all = r_all$cut, f_all = r_all$f,
      n_complete = r_comp$n_seq, cut_complete = r_comp$cut, f_complete = r_comp$f
    )
  }
}

res <- rbindlist(out)
cat("\n\n=== FULL SUMMARY ===\n")
print(as.data.frame(res), row.names = FALSE)

cat("\n=== Pinned-at-1.0 (cut>=0.999) count by iddef x completeness ===\n")
print(res[, .(
  pinned_all = sum(cut_all >= 0.999, na.rm = TRUE),
  pinned_complete = sum(cut_complete >= 0.999, na.rm = TRUE),
  median_f_all = round(median(f_all, na.rm = TRUE), 4),
  median_f_complete = round(median(f_complete, na.rm = TRUE), 4)
), by = iddef])

dir.create("tests/CBSITS_eval/results", showWarnings = FALSE, recursive = TRUE)
fwrite(res, "tests/CBSITS_eval/results/iddef_test.tsv", sep = "\t")
cat("\nWrote tests/CBSITS_eval/results/iddef_test.tsv\n")
