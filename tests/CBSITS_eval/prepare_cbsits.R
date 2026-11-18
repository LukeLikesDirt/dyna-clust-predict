#!/usr/bin/env Rscript
# prepare_cbsits.R — Build a dyna-clust-predict-compatible classification
# file for the CBS-KNAW strain ITS barcode dataset (Vu et al. 2019, Studies
# in Mycology; data/CBSITS/1-s2.0-S0166061618300277-main.pdf).
#
# Inputs:
#   data/CBSITS/CBSITS_classification.fasta  — full-ITS seqs, taxonomy in header
#                                               (k__;p__;c__;o__;f__;g__;s__)
#   data/CBSITS/CBSITS2.fasta                 — ITSx-extracted ITS2-only seqs
#                                               (header = accession only)
#
# Output (mirrors data/ITS2/eukaryome_ITS2.* naming):
#   data/CBSITS_ITS2/cbsits_ITS2.fasta
#   data/CBSITS_ITS2/cbsits_ITS2.classification
#
# Only IDs with a successful ITS2 extraction (present in CBSITS2.fasta) are
# kept, since predict.R requires FASTA and classification IDs to match.
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
})

fasta_its2_in   <- "data/CBSITS/CBSITS2.fasta"
fasta_class_in  <- "data/CBSITS/CBSITS_classification.fasta"
fasta_out       <- "data/CBSITS_ITS2/cbsits_ITS2.fasta"
classification_out <- "data/CBSITS_ITS2/cbsits_ITS2.classification"

dir.create(dirname(fasta_out), recursive = TRUE, showWarnings = FALSE)

# ── ITS2 sequence IDs (successfully extracted by ITSx) ────────────────────
its2_lines  <- readLines(fasta_its2_in)
its2_h_idx  <- which(startsWith(its2_lines, ">"))
its2_ids    <- sub("^>", "", its2_lines[its2_h_idx])
cat("ITS2 sequences available:", length(its2_ids), "\n")

# ── Taxonomy headers (full-ITS file; header carries the CBS taxonomy string) ──
cls_lines <- readLines(fasta_class_in)
cls_h_idx <- which(startsWith(cls_lines, ">"))
cls_headers <- cls_lines[cls_h_idx]

tax_df <- tibble(header = cls_headers) %>%
  mutate(
    id      = str_extract(header, "(?<=^>)[^ ]+"),
    kingdom = str_remove(str_extract(header, "k__([^;]+)"), "k__"),
    phylum  = str_remove(str_extract(header, "p__([^;]+)"), "p__"),
    class   = str_remove(str_extract(header, "c__([^;]+)"), "c__"),
    order   = str_remove(str_extract(header, "o__([^;]+)"), "o__"),
    family  = str_remove(str_extract(header, "f__([^;]+)"), "f__"),
    genus   = str_remove(str_extract(header, "g__([^;]+)"), "g__"),
    sp_raw  = str_remove(str_extract(header, "s__([^ ]+)$"), "s__")
  ) %>%
  mutate(
    # sp_raw is "Genus_epithet[_...]"; strip the leading "Genus_" prefix
    # (matched against this row's own genus field) to recover the epithet,
    # then rebuild "Genus epithet" the same way EUKARYOME's classification
    # files format species (space-separated binomial).
    epithet = if_else(
      !is.na(genus) & str_starts(sp_raw, paste0(genus, "_")),
      str_remove(sp_raw, paste0("^", genus, "_")),
      sp_raw
    ),
    epithet = str_replace_all(epithet, "_", " "),
    species = case_when(
      is.na(genus) | genus %in% c("", "unidentified") ~ "unidentified",
      is.na(epithet) | epithet %in% c("", "unidentified") ~ "unidentified",
      TRUE ~ paste(genus, epithet)
    )
  ) %>%
  select(id, kingdom, phylum, class, order, family, genus, species)

cat("Taxonomy rows parsed:", nrow(tax_df), "\n")

# ── Restrict to IDs with a successful ITS2 extraction ─────────────────────
tax_df <- tax_df %>% filter(id %in% its2_ids)
cat("Rows after restricting to ITS2-extracted IDs:", nrow(tax_df), "\n")

missing_tax <- setdiff(its2_ids, tax_df$id)
if (length(missing_tax) > 0) {
  cat("ITS2 IDs with no taxonomy match (dropped):", length(missing_tax), "\n")
}

if (any(duplicated(tax_df$id))) stop("Duplicate IDs in CBS taxonomy after filtering.")

# ── Write classification ───────────────────────────────────────────────────
write_tsv(tax_df, classification_out)
cat("Classification written to:", classification_out, "\n")

# ── Write matching FASTA (ITS2 seqs restricted to IDs with taxonomy) ──────
keep_ids <- tax_df$id
out_lines <- character(0)
i <- 1
n_written <- 0
while (i <= length(its2_h_idx)) {
  h   <- its2_h_idx[i]
  end <- if (i < length(its2_h_idx)) its2_h_idx[i + 1] - 1L else length(its2_lines)
  id  <- sub("^>", "", its2_lines[h])
  if (id %in% keep_ids) {
    out_lines <- c(out_lines, its2_lines[h:end])
    n_written <- n_written + 1
  }
  i <- i + 1
}
writeLines(out_lines, fasta_out)
cat("FASTA written to:", fasta_out, "(", n_written, "sequences )\n")

# ── Summary ─────────────────────────────────────────────────────────────────
cat("\nTaxonomic composition (identified counts):\n")
source("R/utils.R")
for (rank in c("phylum", "class", "order", "family", "genus", "species")) {
  n_id <- sum(is_identified(tax_df[[rank]]))
  cat(sprintf("  %-8s identified: %-6d unidentified: %d\n",
              rank, n_id, nrow(tax_df) - n_id))
}
