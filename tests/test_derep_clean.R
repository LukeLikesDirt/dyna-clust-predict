################################################################################
# ITS Reference Dataset QC and Dereplication
#
# Steps:
# 1. Load taxonomy and filter to well-represented species
# 2. Inspect sequence length distribution
# 3. Remove sequences outside expected ITS size range
# 4. Remove sequences containing ambiguous nucleotides
# 5. Dereplicate identical sequences
# 6. Resolve taxonomy using LCA
# 7. Write dereplicated fasta and classification files
################################################################################

# ------------------------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------------------------

library(data.table)
library(Biostrings)
library(ggplot2)
library(dplyr)
source("R/utils.R")


# ------------------------------------------------------------------------------
# 1. Read taxonomy
# ------------------------------------------------------------------------------

cls <- fread("data/full_ITS/eukaryome_ITS.classification")


# ------------------------------------------------------------------------------
# 2. Read ITS sequences
# ------------------------------------------------------------------------------

fasta_file <- "data/full_ITS/eukaryome_ITS.fasta"
seqs <- readDNAStringSet(fasta_file)


# ------------------------------------------------------------------------------
# 3. Explore sequence length distribution
# ------------------------------------------------------------------------------

seq_lengths <- width(seqs)

length_df <- data.frame(length = seq_lengths)

summary(seq_lengths)

ggplot(length_df, aes(x = length)) +
  geom_histogram(binwidth = 25, boundary = 0) +
  labs(
    title = "ITS Sequence Length Distribution",
    x = "Sequence length (bp)",
    y = "Count"
  ) +
  xlim(250,1500) +
  theme_minimal()


# ------------------------------------------------------------------------------
# 4. Filter sequences by length
# ------------------------------------------------------------------------------

length_filter <- seq_lengths >= 250 & seq_lengths <= 1500
seqs_length_filtered <- seqs[length_filter]

cat("Original sequences:", length(seqs), "\n")
cat("After length filter:", length(seqs_length_filtered), "\n")


# ------------------------------------------------------------------------------
# 5. Remove sequences with ambiguous nucleotides
# ------------------------------------------------------------------------------

ambiguous_counts <- rowSums(
  letterFrequency(
    seqs_length_filtered,
    letters = c("N","R","Y","S","W","K","M","B","D","H","V")
  )
)

ambiguous_filter <- ambiguous_counts == 0

seqs_clean <- seqs_length_filtered[ambiguous_filter]

cat("Sequences with ambiguous bases:", sum(!ambiguous_filter), "\n")
cat("Sequences after ambiguity filter:", length(seqs_clean), "\n")


# ------------------------------------------------------------------------------
# 6. Dereplicate exact matching sequences
# ------------------------------------------------------------------------------

seq_strings <- as.character(seqs_clean)

derep_table <- data.frame(
  id = names(seqs_clean),
  sequence = seq_strings,
  stringsAsFactors = FALSE
)

# Attach taxonomy
derep_tax <- derep_table %>%
  left_join(cls, by = "id")


# ------------------------------------------------------------------------------
# 7. Resolve taxonomy using LCA
# ------------------------------------------------------------------------------

taxonomy_ranks <- c(
  "kingdom",
  "phylum",
  "class",
  "order",
  "family",
  "genus",
  "species"
)

resolve_lca <- function(df){

  resolved <- list()

  for(rank in taxonomy_ranks){

    vals <- unique(df[[rank]])
    vals <- vals[!vals %in% c("unidentified", NA)]

    if(length(vals) == 1){
      resolved[[rank]] <- vals
    } else {
      resolved[[rank]] <- NA
      break
    }
  }

  as.data.frame(resolved)
}


# Compute LCA taxonomy per unique sequence
resolved_tax <- derep_tax %>%
  group_by(sequence) %>%
  group_modify(~resolve_lca(.x)) %>%
  ungroup()

nrow(resolved_tax)
nrow(derep_tax)
length(seqs_clean)
# ------------------------------------------------------------------------------
# 8. Choose representative sequence per cluster
# ------------------------------------------------------------------------------

# safer alternative to slice()
rep_ids <- derep_table %>%
  group_by(sequence) %>%
  summarise(id = dplyr::first(id), .groups = "drop")


# Attach representative ids
resolved_tax <- resolved_tax %>%
  left_join(rep_ids, by = "sequence")


# ------------------------------------------------------------------------------
# 9. Build dereplicated FASTA
# ------------------------------------------------------------------------------

rep_seqs <- seqs_clean[rep_ids$id]

names(rep_seqs) <- rep_ids$id


# ------------------------------------------------------------------------------
# 10. Build classification table
# ------------------------------------------------------------------------------

classification_out <- resolved_tax %>%
  select(id, all_of(taxonomy_ranks)) %>%
  arrange(id)
# ------------------------------------------------------------------------------
# 11. Write output files
# ------------------------------------------------------------------------------

writeXStringSet(
  rep_seqs,
  "data/full_ITS/eukaryome_ITS_dereplicated.fasta"
)

fwrite(
  classification_out,
  "data/full_ITS/eukaryome_ITS_dereplicated.classification",
  sep = "\t"
)


# ------------------------------------------------------------------------------
# 12. Summary statistics
# ------------------------------------------------------------------------------

cat("\nDereplication summary\n")
cat("----------------------\n")

cat("Sequences after filtering:", length(seqs_clean), "\n")
cat("Unique sequences:", length(rep_seqs), "\n")
cat("Removed duplicates:", length(seqs_clean) - length(rep_seqs), "\n")