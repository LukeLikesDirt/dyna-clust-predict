#!/usr/bin/env Rscript
# append_completeness.R — Append an its_complete column to one or more
# classification files, based on ITSx's raw (pre-dereplication) ITS1 and ITS2
# extraction outputs.
#
# its_complete = TRUE iff the ID's ITS1 AND ITS2 were both extracted by ITSx
# AND both extracted subregions are >= --min_length bp.
#
# "Present in both outputs" alone is not sufficient: it only confirms ITSx
# assigned some label, not that the label represents a genuine full-length
# region. A record with a genuinely short ITS1/ITS2 (min_length protects
# these) is distinct from an ITSx stub call (e.g. a 51bp "ITS1" from a
# severely truncated input record) -- see the dyna-clust-predict cut-off
# investigation for the empirical basis of this distinction.
#
# Completeness is computed once from the RAW ITSx outputs (before
# dereplicate_lca.R collapses identical sequences per subregion), because
# dereplicate_lca()'s classification rebuild only keeps the seven taxonomy
# columns and would silently drop any column appended before it runs. Since
# dereplicate_lca.R groups by exact sequence identity, every ID collapsed into
# one representative shares the same subregion length, so joining
# completeness onto a post-dereplication classification file by ID is safe --
# each surviving representative ID carries its own genuine length.
#
# Usage:
#   Rscript R/append_completeness.R \
#     --its1_raw_fasta      ./tmp/eukaryome_ITS.ITS1.fasta \
#     --its2_raw_fasta      ./tmp/eukaryome_ITS.ITS2.fasta \
#     --min_length          50 \
#     --classification_files data/full_ITS/eukaryome_ITS.classification,data/ITS1/eukaryome_ITS1.classification,data/ITS2/eukaryome_ITS2.classification
#
# Note: This script must be run from the project root directory.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--its1_raw_fasta",
              type = "character", metavar = "FILE",
              help = "Raw ITSx ITS1 output FASTA (pre-dereplication) [required]"),
  make_option("--its2_raw_fasta",
              type = "character", metavar = "FILE",
              help = "Raw ITSx ITS2 output FASTA (pre-dereplication) [required]"),
  make_option("--min_length",
              type = "integer", default = 50L, metavar = "INT",
              help = "Minimum length (bp) for each subregion [default: %default]"),
  make_option("--classification_files",
              type = "character", metavar = "FILE1,FILE2,...",
              help = "Comma-separated classification files to update in place [required]"),
  make_option("--id_col",
              type = "character", default = "id", metavar = "STR",
              help = "ID column name in the classification files [default: %default]")
)

opt <- parse_args(
  OptionParser(
    option_list = option_list,
    usage = "%prog --its1_raw_fasta f1.fasta --its2_raw_fasta f2.fasta --classification_files a.tsv,b.tsv",
    description = "Append an its_complete boolean column to classification files."
  )
)

if (is.null(opt$its1_raw_fasta))       stop("--its1_raw_fasta is required.")
if (is.null(opt$its2_raw_fasta))       stop("--its2_raw_fasta is required.")
if (is.null(opt$classification_files)) stop("--classification_files is required.")
if (!file.exists(opt$its1_raw_fasta))  stop("File not found: ", opt$its1_raw_fasta)
if (!file.exists(opt$its2_raw_fasta))  stop("File not found: ", opt$its2_raw_fasta)

min_length <- opt$min_length
id_col     <- opt$id_col
class_files <- trimws(strsplit(opt$classification_files, ",")[[1]])

# ── Read raw ITSx FASTA lengths (id -> length) ────────────────────────────────

fasta_lengths <- function(path) {
  lines  <- readLines(path)
  h_idx  <- which(startsWith(lines, ">"))
  if (length(h_idx) == 0) return(data.table(id = character(0), len = integer(0)))
  ids <- sub("^>([^ ]+).*", "\\1", lines[h_idx])
  ends <- c(h_idx[-1] - 1L, length(lines))
  lens <- vapply(seq_along(h_idx), function(i) {
    sum(nchar(lines[(h_idx[i] + 1L):ends[i]]))
  }, integer(1))
  data.table(id = ids, len = lens)
}

cat("Reading ITS1 raw lengths from:", opt$its1_raw_fasta, "\n")
its1_len <- fasta_lengths(opt$its1_raw_fasta)
cat("  ", nrow(its1_len), "sequences\n")

cat("Reading ITS2 raw lengths from:", opt$its2_raw_fasta, "\n")
its2_len <- fasta_lengths(opt$its2_raw_fasta)
cat("  ", nrow(its2_len), "sequences\n")

setnames(its1_len, "len", "its1_len")
setnames(its2_len, "len", "its2_len")
merged <- merge(its1_len, its2_len, by = "id")
complete_ids <- merged[its1_len >= min_length & its2_len >= min_length, id]

cat(sprintf(
  "\nComplete (present in both, both >= %dbp): %d of %d present-in-both\n",
  min_length, length(complete_ids), nrow(merged)
))

# ── Update each classification file ───────────────────────────────────────────

for (f in class_files) {
  if (!file.exists(f)) {
    cat("WARNING: classification file not found, skipping:", f, "\n")
    next
  }
  cat("\nUpdating:", f, "\n")
  cls <- fread(f, sep = "\t", header = TRUE, colClasses = "character",
               check.names = FALSE)
  if (!id_col %in% names(cls)) {
    stop("ID column '", id_col, "' not found in ", f)
  }
  cls[, its_complete := get(id_col) %in% complete_ids]
  n_complete <- sum(cls$its_complete)
  cat(sprintf("  %d / %d rows marked its_complete\n", n_complete, nrow(cls)))
  fwrite(cls, f, sep = "\t")
}

cat("\nDone.\n")
