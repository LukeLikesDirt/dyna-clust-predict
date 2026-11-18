#!/usr/bin/env Rscript
# revalidate_hub_filter.R — Re-run the Phase 1c complex-prevalence and
# hub-record diagnosis under the FINALIZED pipeline configuration: 50bp
# completeness floor (2a's actual threshold, not the 100bp used in the first
# pass) and the default --iddef 2 (2c's decision: no change to the default).
#
# Reuses the .sim files already computed in tmp/complex_test/<genus>.sim
# (those were computed on the "present in both, no length floor" superset with
# default iddef=2 -- pairwise scores don't depend on what else is in the
# dataset, so subsetting to the 50bp-floor IDs from that same matrix is valid
# and avoids recomputing vsearch).
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
})

source("R/utils.R")

genera <- c("Fusarium", "Alternaria", "Trichoderma", "Cladosporium",
           "Penicillium", "Aspergillus")

refined_ids <- fread("/tmp/refined_complete_ids_50bp.txt", header = FALSE)$V1

cls <- fread("data/full_ITS/eukaryome_ITS.classification")

hub_threshold <- 3L  # max distinct partner species before a record is flagged as a hub

out <- list()
for (g in genera) {
  sim_file <- file.path("tmp/complex_test", paste0(g, ".sim"))
  if (!file.exists(sim_file)) { cat("missing sim for", g, "\n"); next }

  sub <- cls[genus == g & is_identified(species) & id %in% refined_ids]
  if (nrow(sub) < 30) { cat("skip (too few 50bp-floor seqs):", g, "\n"); next }

  sim <- fread(sim_file, header = FALSE, col.names = c("q", "t", "pid"))
  sim <- sim[q %in% sub$id & t %in% sub$id]

  m <- merge(sim, sub[, .(id, species)], by.x = "q", by.y = "id"); setnames(m, "species", "sq")
  m <- merge(m, sub[, .(id, species)], by.x = "t", by.y = "id"); setnames(m, "species", "st")

  seqcount <- sub[, .(N = .N), by = species]
  spp <- sort(unique(sub$species))

  # ── Baseline: unguarded single-linkage complex detection ──────────────────
  e_full <- unique(m[pid >= 100 & sq != st, .(a = pmin(sq, st), b = pmax(sq, st))])
  g_full <- graph_from_data_frame(e_full, directed = FALSE, vertices = data.frame(name = spp))
  cm_full <- components(g_full)
  dt_full <- data.table(species = spp, comp = cm_full$membership[spp])
  dt_full <- merge(dt_full, seqcount, by = "species", all.x = TRUE); dt_full[is.na(N), N := 0L]
  setorder(dt_full, comp, -N)
  dt_full[, keep := seq_len(.N) == 1L, by = comp]
  csz_full <- dt_full[, .N, by = comp]
  n_complexes_full   <- nrow(csz_full[N > 1])
  largest_full       <- if (nrow(csz_full) > 0) max(csz_full$N) else 0L
  species_removed_full <- sum(!dt_full$keep)

  # ── Hub-guarded: exclude records whose 100%-identity partners span too many species ──
  e_seq <- m[pid >= 100 & sq != st]
  hub_score <- rbindlist(list(
    e_seq[, .(id = q, partner_species = st)],
    e_seq[, .(id = t, partner_species = sq)]
  ))[, .(n_distinct_partner_species = uniqueN(partner_species)), by = id]
  bad_ids <- hub_score[n_distinct_partner_species >= hub_threshold, id]

  m2 <- m[!(q %in% bad_ids) & !(t %in% bad_ids)]
  e_guard <- unique(m2[pid >= 100 & sq != st, .(a = pmin(sq, st), b = pmax(sq, st))])
  g_guard <- graph_from_data_frame(e_guard, directed = FALSE, vertices = data.frame(name = spp))
  cm_guard <- components(g_guard)
  dt_guard <- data.table(species = spp, comp = cm_guard$membership[spp])
  dt_guard <- merge(dt_guard, seqcount, by = "species", all.x = TRUE); dt_guard[is.na(N), N := 0L]
  setorder(dt_guard, comp, -N)
  dt_guard[, keep := seq_len(.N) == 1L, by = comp]
  csz_guard <- dt_guard[, .N, by = comp]
  n_complexes_guard    <- nrow(csz_guard[N > 1])
  largest_guard        <- if (nrow(csz_guard) > 0) max(csz_guard$N) else 0L
  species_removed_guard <- sum(!dt_guard$keep)

  cat(sprintf(
    "%-13s n_seq=%-5d n_species=%-4d hub_records=%-3d | UNGUARDED: complexes=%-3d largest=%-3d removed=%-3d (%.1f%%) | HUB-FILTERED: complexes=%-3d largest=%-3d removed=%-3d (%.1f%%)\n",
    g, nrow(sub), length(spp), length(bad_ids),
    n_complexes_full, largest_full, species_removed_full, round(100*species_removed_full/length(spp),1),
    n_complexes_guard, largest_guard, species_removed_guard, round(100*species_removed_guard/length(spp),1)
  ))

  out[[g]] <- data.table(
    genus = g, n_seq = nrow(sub), n_species = length(spp), hub_records = length(bad_ids),
    complexes_unguarded = n_complexes_full, largest_unguarded = largest_full,
    removed_unguarded = species_removed_full,
    complexes_guarded = n_complexes_guard, largest_guarded = largest_guard,
    removed_guarded = species_removed_guard
  )
}

res <- rbindlist(out)
cat("\n=== SUMMARY: 50bp completeness floor, default iddef=2 ===\n")
print(as.data.frame(res), row.names = FALSE)

dir.create("tests/CBSITS_eval/results", showWarnings = FALSE, recursive = TRUE)
fwrite(res, "tests/CBSITS_eval/results/hub_filter_revalidation_50bp.tsv", sep = "\t")
cat("\nWrote tests/CBSITS_eval/results/hub_filter_revalidation_50bp.tsv\n")
