
# Required packages
require(data.table)
require(dplyr)

# Source utility functions
source("R/utils.R")

# Read in the classification data
taxa <- fread("data/full_ITS/eukaryome_ITS.classification")

# Identified species
species <- taxa %>%
  filter(
    is_identified(species)
    )

# Count the number of species per kingdom
species %>%
  group_by(kingdom) %>%
  summarise(n_species = n_distinct(species)) %>%
  print(n = Inf)

# Count the number of species per phylum
species %>%
  group_by(kingdom, phylum) %>%
  summarise(n_species = n_distinct(species)) %>%
  print(n = Inf)

# Count the number of species per order
species %>%
  group_by(kingdom, phylum, order) %>%
  summarise(n_species = n_distinct(species)) %>%
  print(n = 600)

# Count the number of sequences per species
species %>%
  group_by(kingdom, phylum, species) %>%
  summarise(
    n_secs = n()
    ) %>%
  filter(n_secs == 1)

species %>%filter ( kingdom=="_mitochondrion")
