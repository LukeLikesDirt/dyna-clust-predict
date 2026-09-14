#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=0-01:00:00
#SBATCH --partition=short
#SBATCH --output=logs/%x.%j.out

# Script name:  05_prepare_subsets.sh
# Description:  Prepare prediction ID files for all three ITS regions (full ITS,
#               ITS1, ITS2) by running subset.R once per region. Produces one
#               ID file per unique-sequence rank (STEP 1) and one ID file per
#               valid (target x parent) rank combination (STEP 2).
#
#               Reads each region's *_nocomplex.classification (written by
#               04_remove_complexes.sh), not the pre-complex-removal file --
#               run 04 before this script.
# Note:         This script must be run from the project root directory.

# =============================================================================
# PARAMETER SETUP
# =============================================================================

readonly SUBSET="./R/subset.R"

# Filter constants
readonly MIN_SUBGROUPS=10
readonly MIN_SEQUENCES=30
readonly MAX_SEQUENCES=20000
readonly MAX_PROPORTION=0.5
readonly MAX_KINGDOM_PROPORTION=0.5
# Caps how many sequences any single SPECIES may contribute before the
# max_sequences downsample (a no-op at every coarser target rank -- see
# species_only() in R/subset.R). A species with many sequences and real
# internal diversity can't all be pairwise-identical near threshold 1.0, so
# it creates useful downward pressure on the optimal cutoff -- but only if
# enough of its sequences survive to reveal that diversity, and a small
# number of species are oversampled in EUKARYOME for non-biological reasons
# (e.g. medically/agriculturally important fungi): 3.3% of Fungi species
# carry 26% of all sequences. 10 keeps meaningful within-species signal
# while bounding how much any one such species can dominate the
# sequence-count-weighted F-measure.
readonly MAX_SEQS_PER_GROUP=10
# Caps the fraction of SPECIES (by count) that may be singletons (exactly 1
# sequence) -- also species-target only. A singleton scores a free Dice =
# 1.0 at threshold 1.0, pinning singleton-dominated datasets' cutoffs near
# 1.0 by construction regardless of biology. 0.7 was chosen empirically on
# kingdom Fungi (tests/CBSITS_eval/test_fungi_grid.R): combined with
# MAX_SEQS_PER_GROUP=10, it lands the cutoff near 0.97, a conventional and
# defensible value for fungal ITS, without discarding as much real data as
# the tighter 0.5 cap tested earlier. Weaker self-predictions this pair
# pushes toward 1.0 (e.g. a well-sampled genus that no longer needs
# downsampling to demonstrate its own diversity) are expected to be rescued
# by consolidate_cutoffs.R's confidence-ranked fallback to a higher-rank
# ancestor, which carries far more multi-sequence evidence at production
# scale and clears that guard easily.
readonly MAX_SINGLETON_PROPORTION=0.7

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

echo "Activating conda environment..."
source ~/.bashrc
conda activate dyna_clust_predict

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if [[ ! -f "$SUBSET" ]]; then
    echo "ERROR: R script not found: $SUBSET" >&2
    exit 1
fi

# =============================================================================
# HELPER FUNCTION
# =============================================================================

# run_subset <classification_in> <output_dir> <region_label>
run_subset() {
  local classification_in="$1"
  local output_dir="$2"
  local label="$3"

  if [[ ! -f "$classification_in" ]]; then
    echo "WARNING: Classification not found for region '$label', skipping: $classification_in" >&2
    return 0
  fi

  echo ""
  echo "--- Region: $label ---"
  echo "Classification : $classification_in"
  echo "Output dir     : $output_dir"
  echo "$(date)"

  Rscript "$SUBSET" \
    --classification_in "$classification_in" \
    --output_dir        "$output_dir" \
    --min_subgroups     "$MIN_SUBGROUPS" \
    --min_sequences     "$MIN_SEQUENCES" \
    --max_sequences     "$MAX_SEQUENCES" \
    --max_proportion    "$MAX_PROPORTION" \
    --max_kingdom_proportion "$MAX_KINGDOM_PROPORTION" \
    --max_seqs_per_group "$MAX_SEQS_PER_GROUP" \
    --max_singleton_proportion "$MAX_SINGLETON_PROPORTION"

  if [[ $? -ne 0 ]]; then
    echo "ERROR: subset.R failed for region '$label'." >&2
    return 1
  fi
  echo "Finished region '$label' at: $(date)"
}

# =============================================================================
# PREPARE SUBSETS FOR ALL THREE REGIONS
# =============================================================================

echo ""
echo "=== PREPARING SUBSETS ==="
echo "min_subgroups  : $MIN_SUBGROUPS"
echo "min_sequences  : $MIN_SEQUENCES"
echo "max_sequences  : $MAX_SEQUENCES"
echo "max_proportion : $MAX_PROPORTION"
echo "max_kingdom_proportion : $MAX_KINGDOM_PROPORTION"
echo "max_seqs_per_group : $MAX_SEQS_PER_GROUP"
echo "max_singleton_proportion : $MAX_SINGLETON_PROPORTION"

# 1. Full ITS
run_subset \
    "./data/full_ITS/eukaryome_ITS_nocomplex.classification" \
    "./data/full_ITS" \
    "full_ITS"

# 2. ITS1
run_subset \
    "./data/ITS1/eukaryome_ITS1_nocomplex.classification" \
    "./data/ITS1" \
    "ITS1"

# 3. ITS2
run_subset \
    "./data/ITS2/eukaryome_ITS2_nocomplex.classification" \
    "./data/ITS2" \
    "ITS2"

echo ""
echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="
echo "$(date)"
echo ""
echo "ID files written to data/full_ITS, data/ITS1, data/ITS2"
echo "  STEP 1 unique-sequence IDs : <rank>_unique_id.txt"
echo "  STEP 2 prediction IDs      : <target>_pred_id_<parent>.txt"

conda deactivate
