#!/usr/bin/env Rscript
# test_reciprocal_support.R — Compare the hub-exclusion safeguard (currently
# implemented in R/remove_complexes.R) against the untested reciprocal/
# multi-pair-support alternative: require >= --min_pairs independent
# sequence-pairs at the threshold identity before accepting a species-pair
# edge into the complex graph, instead of excluding individual hub sequences.
#
# Rationale for reciprocal support: a single mislabeled/contaminated record
# can produce exactly one spurious 100%-identity pair between two otherwise
# unrelated species. Requiring >=2 independent pairs makes a lone bad record
# insufficient to create a false edge, without needing to identify which
# specific sequence is the outlier.
#
# Runs both mechanisms on the same 6 genera, same 50bp-completeness-floor
# input, same cached .sim files (default iddef=2) already used to validate
# hub-exclusion, so the comparison is apples-to-apples.
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

hub_threshold  <- 3L  # current R/remove_complexes.R default
min_pairs_opts <- c(2L, 3L)  # reciprocal-support thresholds to test

out <- list()
for (g in genera) {
  sim_file <- file.path("tmp/complex_test", paste0(g, ".sim"))
  if (!file.exists(sim_file)) { cat("missing sim for", g, "\n"); next }

  sub <- cls[genus == g & is_identified(species) & id %in% refined_ids]
  if (nrow(sub) < 30) { cat("skip:", g, "\n"); next }

  sim <- fread(sim_file, header = FALSE, col.names = c("q", "t", "pid"))
  sim <- sim[q %in% sub$id & t %in% sub$id]

  m <- merge(sim, sub[, .(id, species)], by.x = "q", by.y = "id"); setnames(m, "species", "sq")
  m <- merge(m, sub[, .(id, species)], by.x = "t", by.y = "id"); setnames(m, "species", "st")
  m <- m[pid >= 100 & sq != st]

  seqcount <- sub[, .(N = .N), by = species]
  spp <- sort(unique(sub$species))

  build_and_resolve <- function(edges) {
    g_graph <- graph_from_data_frame(edges, directed = FALSE, vertices = data.frame(name = spp))
    cm <- components(g_graph)
    dt <- data.table(species = spp, comp = cm$membership[spp])
    dt <- merge(dt, seqcount, by = "species", all.x = TRUE); dt[is.na(N), N := 0L]
    setorder(dt, comp, -N, species)
    dt[, keep := seq_len(.N) == 1L, by = comp]
    csz <- dt[, .N, by = comp]
    list(n_complexes = nrow(csz[N > 1]),
        largest = if (nrow(csz) > 0) max(csz$N) else 0L,
        removed = sum(!dt$keep))
  }

  # ── Method 1: hub-exclusion (current implementation) ───────────────────────
  hub_score <- rbindlist(list(
    m[, .(id = q, partner_species = st)],
    m[, .(id = t, partner_species = sq)]
  ))[, .(n_distinct_partner_species = uniqueN(partner_species)), by = id]
  hub_ids <- hub_score[n_distinct_partner_species >= hub_threshold, id]
  m_hub <- m[!(q %in% hub_ids) & !(t %in% hub_ids)]
  e_hub <- unique(m_hub[, .(a = pmin(sq, st), b = pmax(sq, st))])
  r_hub <- build_and_resolve(e_hub)

  # ── Method 2: reciprocal/multi-pair support, at each min_pairs threshold ───
  pair_counts <- unique(m[, .(q, t, sq, st)])[, .(a = pmin(sq, st), b = pmax(sq, st)), by = .(q, t)]
  species_pair_support <- pair_counts[, .N, by = .(a, b)]

  r_recip <- lapply(min_pairs_opts, function(mp) {
    e_recip <- species_pair_support[N >= mp, .(a, b)]
    build_and_resolve(e_recip)
  })
  names(r_recip) <- paste0("recip", min_pairs_opts)

  cat(sprintf(
    "%-13s HUB-EXCL(>=%d): complexes=%-3d largest=%-3d removed=%-3d | RECIP(>=2 pairs): complexes=%-3d largest=%-3d removed=%-3d | RECIP(>=3 pairs): complexes=%-3d largest=%-3d removed=%-3d\n",
    g, hub_threshold, r_hub$n_complexes, r_hub$largest, r_hub$removed,
    r_recip$recip2$n_complexes, r_recip$recip2$largest, r_recip$recip2$removed,
    r_recip$recip3$n_complexes, r_recip$recip3$largest, r_recip$recip3$removed
  ))

  out[[g]] <- data.table(
    genus = g, n_species = length(spp),
    hub_complexes = r_hub$n_complexes, hub_largest = r_hub$largest, hub_removed = r_hub$removed,
    recip2_complexes = r_recip$recip2$n_complexes, recip2_largest = r_recip$recip2$largest, recip2_removed = r_recip$recip2$removed,
    recip3_complexes = r_recip$recip3$n_complexes, recip3_largest = r_recip$recip3$largest, recip3_removed = r_recip$recip3$removed
  )
}

res <- rbindlist(out)
cat("\n=== SUMMARY ===\n")
print(as.data.frame(res), row.names = FALSE)

dir.create("tests/CBSITS_eval/results", showWarnings = FALSE, recursive = TRUE)
fwrite(res, "tests/CBSITS_eval/results/reciprocal_support_comparison.tsv", sep = "\t")
cat("\nWrote tests/CBSITS_eval/results/reciprocal_support_comparison.tsv\n")
