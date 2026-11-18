################################################################################
# ITS Reference Dataset QC and Dereplication
#
# Steps:
# 1. Load taxonomy and filter to well-represented species
# 2. Inspect sequence length distribution
# 3. Remove sequences outside expected ITS size range
# 4. Remove sequences containing ambiguous nucleotides
# 5. Dereplicate identical sequences
# 6. Identify identical sequences with:
#    a) same taxonomy
#    b) conflicting taxonomy
################################################################################


# ------------------------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(Biostrings)
library(ggplot2)
source("R/utils.R")

# ------------------------------------------------------------------------------
# 1. Read taxonomy and filter species with sufficient representation
# ------------------------------------------------------------------------------
cls <- fread("data/full_ITS/eukaryome_ITS.classification")

cls_filtered <- cls %>%
  filter(
    is_identified(species),                # keep only identified species
    is_identified(kingdom)
  )

print(cls_filtered)

# ------------------------------------------------------------------------------
# Summary of dataset by kingdom
# ------------------------------------------------------------------------------

cls_filtered %>%
  group_by(kingdom, phylum) %>%
  summarise(
    n_species = n_distinct(species),
    n_seqs = n()
  ) %>%
  group_by(kingdom, phylum) %>%
  arrange(kingdom, phylum, n_species) %>%
  print(n = Inf)

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
  xlim(250, 1500) +
  theme_minimal()

# ------------------------------------------------------------------------------
# 4. Filter sequences by length
# ------------------------------------------------------------------------------

length_filter <- seq_lengths >= 250 & seq_lengths <= 1500

seqs_length_filtered <- seqs[length_filter]

cat("Original sequences:", length(seqs), "\n")
cat("After length filter:", length(seqs_length_filtered), "\n")

# ------------------------------------------------------------------------------
# 5. Identify sequences containing ambiguous nucleotides
#
# Ambiguous IUPAC bases include:
# N R Y S W K M B D H V
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

# Convert sequences to character strings for grouping
seq_strings <- as.character(seqs_clean)

derep_table <- data.frame(
  id = names(seqs_clean),
  sequence = seq_strings,
  stringsAsFactors = FALSE
)

# Group identical sequences
derep_groups <- derep_table %>%
  group_by(sequence) %>%
  mutate(cluster_size = n()) %>%
  ungroup()

# ------------------------------------------------------------------------------
# 7. Attach taxonomy information
# ------------------------------------------------------------------------------

derep_tax <- derep_groups %>%
  left_join(cls, by = "id")

# ------------------------------------------------------------------------------
# 8. Identify identical sequences with taxonomy conflicts
# ------------------------------------------------------------------------------

taxonomy_consistency <- derep_tax %>%
  group_by(sequence) %>%
  summarise(
    n_sequences = n(),
    n_species = n_distinct(species),
    species_list = paste(unique(species), collapse = ";"),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# Exact matches with identical taxonomy
# ------------------------------------------------------------------------------

exact_same_tax <- taxonomy_consistency %>%
  filter(n_sequences > 1 & n_species == 1)

# ------------------------------------------------------------------------------
# Exact matches with conflicting taxonomy
# ------------------------------------------------------------------------------

exact_conflicting_tax <- taxonomy_consistency %>%
  filter(n_sequences > 1 & n_species > 1)

# ------------------------------------------------------------------------------
# Resolving duplicates
# ------------------------------------------------------------------------------

# These can be dereplicated to a single representative sequence without losing taxonomic information
exact_same_tax %>%
  select(n_sequences, n_species, species_list) %>%
  filter(species_list != "unidentified")

# Identified species that conflict need to be resolved to the lowest common taxonomic level (e.g. genus or family)
exact_conflicting_tax %>%
  filter(!grepl("unidentified", species_list)) %>%
  select(n_sequences, n_species, species_list) %>%
  print(n = 10)

# For these cases, sequences need to be resolved to the lowest common taxonomic level (e.g. genus or family) to avoid conflicts.
# However, if they one is considered "unidentified" at a given rank (e.g. family), and the other is identified at that rank, then the identified taxonomy can be retained if the first rank where both are identified matches 
# (similar to the "unidentified" case below).
exact_same_tax %>%
  filter(species_list == "unidentified") %>%
  select(n_sequences, n_species, species_list) %>%
  print(n = 10)
exact_conflicting_tax %>%
  filter(grepl("unidentified", species_list)) %>%
  select(n_sequences, n_species, species_list) %>%
  print(n = 10)

# ------------------------------------------------------------------------------
# Summary statistics
# ------------------------------------------------------------------------------

cat("\nDereplication summary\n")
cat("----------------------\n")
cat("Unique sequences:", nrow(taxonomy_consistency), "\n")
cat("Exact matches (same taxonomy):", nrow(exact_same_tax), "\n")
cat("Exact matches (conflicting taxonomy):", nrow(exact_conflicting_tax), "\n")

# ------------------------------------------------------------------------------
# Optional: inspect problematic clusters
# ------------------------------------------------------------------------------

exact_conflicting_tax %>%
  select(n_sequences, n_species, species_list) %>%
  print(n = 100)
