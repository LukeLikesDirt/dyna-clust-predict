#!/usr/bin/env Rscript
# test_complex_prevalence.R — Phase 1c: real complex prevalence once the
# completeness filter removes truncation artifacts, using the DEFAULT
# identity definition (--iddef 2) -- per the Phase 1a/1b finding that
# switching --iddef costs more F-measure than it recovers, while the
# completeness filter alone already un-pins every real fungal genus tested.
#
# The 64.7% species-removal figure measured earlier (Fusarium, iddef=2, ALL
# sequences including partial/truncated ones) was the artifact baseline, not
# the expected real prevalence. This reruns single-linkage complex detection
# at 100% identity restricted to complete-span sequences.
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
})

source("R/utils.R")

tmp_dir <- "tmp/complex_test"
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

genera <- c("Fusarium", "Alternaria", "Trichoderma", "Cladosporium",
           "Penicillium", "Aspergillus")

fasta_all <- "data/full_ITS/eukaryome_ITS.fasta"
class_all <- "data/full_ITS/eukaryome_ITS.classification"

cat("Building complete-span ID set...\n")
get_ids <- function(f) sub("^>", "", grep("^>", readLines(f), value = TRUE))
its1_ids <- get_ids("data/ITS1/eukaryome_ITS1.fasta")
its2_ids <- get_ids("data/ITS2/eukaryome_ITS2.fasta")
complete_ids <- intersect(its1_ids, its2_ids)
cat("Complete-span IDs:", length(complete_ids), "\n\n")

cls <- fread(class_all)

out <- list()
for (g in genera) {
  sub <- cls[genus == g & is_identified(species)]
  sub <- sub[id %in% complete_ids]
  if (nrow(sub) < 30) { cat("skip (too few complete-span seqs):", g, "\n"); next }

  ids_file <- file.path(tmp_dir, paste0(g, ".ids"))
  fa_file  <- file.path(tmp_dir, paste0(g, ".fasta"))
  sim_file <- file.path(tmp_dir, paste0(g, ".sim"))
  writeLines(sub$id, ids_file)
  system(sprintf("seqkit grep -f '%s' '%s' > '%s' 2>/dev/null", ids_file, fasta_all, fa_file))
  system(sprintf(paste("vsearch --allpairs_global '%s' --acceptall --userout '%s'",
                       "--userfields query+target+id --threads 10"),
                 fa_file, sim_file), ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (!file.exists(sim_file) || file.size(sim_file) == 0) { cat("no sim:", g, "\n"); next }

  sim <- fread(sim_file, header = FALSE, col.names = c("q", "t", "pid"))
  seqcount <- sub[, .(N = .N), by = species]
  spp <- sort(unique(sub$species))

  m <- merge(sim, sub[, .(id, species)], by.x = "q", by.y = "id"); setnames(m, "species", "sq")
  m <- merge(m, sub[, .(id, species)], by.x = "t", by.y = "id"); setnames(m, "species", "st")

  e <- unique(m[pid >= 100 & sq != st, .(a = pmin(sq, st), b = pmax(sq, st))])
  g_graph <- graph_from_data_frame(e, directed = FALSE, vertices = data.frame(name = spp))
  cm <- components(g_graph)
  dt <- data.table(species = spp, comp = cm$membership[spp])
  dt <- merge(dt, seqcount, by = "species", all.x = TRUE); dt[is.na(N), N := 0L]
  setorder(dt, comp, -N)
  dt[, keep := seq_len(.N) == 1L, by = comp]
  csz <- dt[, .N, by = comp]

  n_species <- length(spp)
  n_complexes <- nrow(csz[N > 1])
  largest_complex <- max(csz$N)
  species_removed <- sum(!dt$keep)
  seqs_removed <- sum(dt[keep == FALSE]$N)

  cat(sprintf("%-13s complete-span n_seq=%-5d n_species=%-4d complexes=%-3d largest=%-3d species_removed=%-3d (%.1f%%) seqs_removed=%d (%.1f%%)\n",
              g, nrow(sub), n_species, n_complexes, largest_complex,
              species_removed, round(100*species_removed/n_species,1),
              seqs_removed, round(100*seqs_removed/sum(dt$N),1)))

  out[[g]] <- data.table(genus = g, n_seq = nrow(sub), n_species = n_species,
                         n_complexes = n_complexes, largest_complex = largest_complex,
                         species_removed = species_removed,
                         pct_species_removed = round(100*species_removed/n_species,1),
                         seqs_removed = seqs_removed,
                         pct_seqs_removed = round(100*seqs_removed/sum(dt$N),1))

  # save the complex membership for manifest-style review
  dt[, genus := g]
  fwrite(dt[comp %in% csz[N>1]$comp][order(comp,-N)],
         file.path(tmp_dir, paste0(g, "_complexes.tsv")), sep="\t")
}

res <- rbindlist(out)
cat("\n=== SUMMARY: complex prevalence, complete-span sequences, default iddef ===\n")
print(as.data.frame(res), row.names = FALSE)

dir.create("tests/CBSITS_eval/results", showWarnings = FALSE, recursive = TRUE)
fwrite(res, "tests/CBSITS_eval/results/complex_prevalence_test.tsv", sep = "\t")
cat("\nWrote tests/CBSITS_eval/results/complex_prevalence_test.tsv\n")
cat("Per-genus complex membership: tmp/complex_test/<genus>_complexes.tsv\n")
