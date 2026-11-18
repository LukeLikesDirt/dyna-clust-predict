#!/usr/bin/env Rscript
# create_multi_test_datasets.R — Create N independent 2000-seq test datasets
# from the full ITS1 eukaryome data
#
# Usage: Rscript create_multi_test_datasets.R <fasta> <classification> <output_dir> <n_datasets> <n_seqs>

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
input_fasta <- args[1]
input_class <- args[2]
output_dir  <- args[3]
n_datasets  <- as.integer(args[4])
n_seqs      <- as.integer(args[5])

cls <- fread(input_class, sep = "\t", header = TRUE, quote = "",
             fill = TRUE, na.strings = "", check.names = FALSE)
all_ids <- cls$id
cat(sprintf("Total IDs: %d\n", length(all_ids)))

# Sample without replacement across datasets to ensure independence
set.seed(123)
if (n_datasets * n_seqs > length(all_ids)) {
  stop("Not enough sequences for requested datasets")
}

shuffled <- sample(all_ids)

for (i in seq_len(n_datasets)) {
  start_idx <- (i - 1) * n_seqs + 1
  end_idx   <- i * n_seqs
  ds_ids    <- shuffled[start_idx:end_idx]

  ds_dir <- file.path(output_dir, sprintf("dataset_%d", i))
  dir.create(ds_dir, showWarnings = FALSE, recursive = TRUE)

  # Write IDs
  id_file <- file.path(ds_dir, "ids.txt")
  writeLines(ds_ids, id_file)

  # Write classification subset
  cls_sub <- cls[cls$id %in% ds_ids, ]
  cls_file <- file.path(ds_dir, "test.classification")
  fwrite(cls_sub, cls_file, sep = "\t", quote = FALSE)

  cat(sprintf("Dataset %d: %d sequences written to %s\n", i, length(ds_ids), ds_dir))
}

cat("Done.\n")
