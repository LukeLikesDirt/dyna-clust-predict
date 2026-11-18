#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=0-01:00:00
#SBATCH --partition=short
#SBATCH --output=logs/%x.%j.out

# Script name:  05_prepare_subsets.sh
# Description:  Write the sequence ID files that define each prediction
#               dataset, for full ITS, ITS1 and ITS2 (R/subset.R):
#               - STEP 1: <rank>_unique_id.txt, one per rank
#               - STEP 2: <target>_pred_id_<parent>.txt, one per target x
#                 parent rank combination, plus <target>_pred_id_global.txt
#
# Note:         This script must be run from the project root directory,
#               after 04.

# =============================================================================
# PARAMETER SETUP
# =============================================================================

readonly SUBSET="./R/subset.R"

# Minimum number of child taxa per dataset.
readonly MIN_SUBGROUPS=10
# Minimum number of sequences per dataset.
readonly MIN_SEQUENCES=30
# Datasets larger than this are downsampled by randomly selecting whole child
# taxa (e.g. whole species). Must not exceed predict.R's --max_seq_no (20,000).
readonly MAX_SEQUENCES=20000
# Maximum proportion of sequences a child taxon may contribute to the parent
# group.
readonly MAX_PROPORTION=0.7
# Maximum proportion of sequences a kingdom may contribute to global datasets.
readonly MAX_KINGDOM_PROPORTION=0.5
# Maximum sequences per species, applied before the MAX_SEQUENCES downsample.
# Species datasets only.
readonly MAX_SEQS_PER_GROUP=10
# Maximum proportion of singleton species a dataset may contain. A singleton
# scores a free Dice = 1.0 at threshold 1.0, so a singleton-dominated dataset
# has its optimum pinned near 1.0. Species datasets only.
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
