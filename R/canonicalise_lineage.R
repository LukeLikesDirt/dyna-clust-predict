#!/usr/bin/env Rscript
# canonicalise_lineage.R — Make the classification table a tree: every named
# taxon gets exactly one parent.
#
# Run by 01_reformat_ITS.sh straight after reformat.R, so every file built from
# the classification downstream (dereplication, ITSx subregions, complex
# removal, subsets, predictions, consolidation) sees one lineage per name.
# Without it, consolidate_cutoffs.R warns that some names "have an inconsistent
# parent ... across the classification file" and looks their ancestors up by
# the modal parent, and dereplicate_lca.R truncates the LCA of identical
# sequences at the conflicting rank.
#
# The rule (canonicalise_lineage() in R/utils.R): for each named label, the
# modal parent among identified parents, else among incertae sedis parents;
# ties to the alphabetically first. Records under a marker kingdom
# (_mitochondrion, _pseudogene ...) are left untouched. Every label that had
# more than one parent is listed in --manifest_out with the parent chosen, so a
# new EUKARYOME release shows up as a reviewable diff.
#
# Usage:
#   Rscript R/canonicalise_lineage.R \
#     --classification_in  tmp/eukaryome_ITS.classification \
#     --classification_out tmp/eukaryome_ITS.classification \
#     --manifest_out       data/lineage_canonicalisation.txt
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

source("R/utils.R")

# ── Arguments ─────────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--classification_in",
              type = "character", metavar = "FILE",
              help = "Tab-delimited classification (id + kingdom..species) [required]"),
  make_option("--classification_out",
              type = "character", metavar = "FILE",
              help = "Output classification; may be the same path as the input [required]"),
  make_option("--manifest_out",
              type = "character", default = "data/lineage_canonicalisation.txt", metavar = "FILE",
              help = "Labels that had more than one parent, and the parent chosen [default: %default]")
)

opt <- parse_args(
  OptionParser(
    option_list = option_list,
    usage       = "%prog --classification_in in.tsv --classification_out out.tsv",
    description = "Give every named taxon in a classification table exactly one parent."
  )
)

if (is.null(opt$classification_in))  stop("--classification_in is required.")
if (is.null(opt$classification_out)) stop("--classification_out is required.")
if (!file.exists(opt$classification_in)) stop("File not found: ", opt$classification_in)

for (d in unique(dirname(c(opt$classification_out, opt$manifest_out)))) {
  if (!dir.exists(d) && d != ".") dir.create(d, recursive = TRUE)
}

# ── Canonicalise ──────────────────────────────────────────────────────────────

cat("[canonicalise] Loading classification from:", opt$classification_in, "\n")
# All text, no quote parsing and no NA strings: a label is kept exactly as written.
cls <- fread(opt$classification_in, sep = "\t", header = TRUE, quote = "",
             colClasses = "character", na.strings = NULL)
cat(sprintf("[canonicalise] %d records\n", nrow(cls)))

before <- find_nesting_violations(cls[!startsWith(kingdom, "_")])
res    <- canonicalise_lineage(cls)

cat(sprintf("[canonicalise] labels with more than one parent: %d before, 0 after\n", nrow(before)))
if (nrow(res$manifest)) {
  by_rank <- res$manifest[, .N, by = rank][order(match(rank, rank_hierarchy))]
  cat("[canonicalise]   by rank:", paste(sprintf("%s %d", by_rank$rank, by_rank$N), collapse = ", "), "\n")
}
cat(sprintf("[canonicalise] records re-parented: %d (%.3f%%)\n",
            res$n_rows_changed, 100 * res$n_rows_changed / max(1L, nrow(cls))))
not_modal <- res$manifest[n_records_changed > n_records / 2]
cat(sprintf("[canonicalise] labels whose chosen parent is not the majority parent: %d",
            nrow(not_modal)),
    "(an identified parent beats an incertae sedis or missing one)\n")
cat("[canonicalise] named labels under a missing parent with no identified parent anywhere",
    "(left as they are):",
    paste(sprintf("%s %d", names(res$n_orphans), res$n_orphans), collapse = ", "), "\n")

# ── Write outputs ─────────────────────────────────────────────────────────────

fwrite(res$classification, opt$classification_out, sep = "\t", quote = FALSE, na = "")
fwrite(res$manifest[order(match(rank, rank_hierarchy), child)], opt$manifest_out,
       sep = "\t", quote = FALSE, na = "")

cat("[canonicalise] Classification written to:", opt$classification_out, "\n")
cat("[canonicalise] Manifest written to:      ", opt$manifest_out, "\n")
