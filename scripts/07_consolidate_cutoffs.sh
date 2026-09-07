#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0-01:00:00
#SBATCH --partition=short
#SBATCH --output=logs/%x.%j.out

# Script name:  07_consolidate_cutoffs.sh
# Description:  Fill gaps and repair monotonicity in each region's nested
#               cutoff table via R/consolidate_cutoffs.R. Pure table work (no
#               vsearch, no sequence I/O), so this is fast and single-threaded
#               in effect -- the SBATCH resources above are generous headroom,
#               not a requirement.
#
#               Requires 06a_predict_cutoffs.sh / 06c_predict_cutoffs_parallel.sh
#               to have already produced eukaryome.cutoffs.json.txt for the
#               region, INCLUDING the global (no --higher_rank) rows -- this
#               script's fallback chain has no top-level anchor without them.
#
# Note:         This script must be run from the project root directory.

# =============================================================================
# PARAMETER SETUP
# =============================================================================

readonly CONSOLIDATE="./R/consolidate_cutoffs.R"
readonly PREFIX="eukaryome"

REGION_LABELS=("full_ITS" "ITS1" "ITS2")

declare -A REGION_CUTOFFS_IN
REGION_CUTOFFS_IN["full_ITS"]="./data/full_ITS/${PREFIX}.cutoffs.json.txt"
REGION_CUTOFFS_IN["ITS1"]="./data/ITS1/${PREFIX}.cutoffs.json.txt"
REGION_CUTOFFS_IN["ITS2"]="./data/ITS2/${PREFIX}.cutoffs.json.txt"

declare -A REGION_CLASS
REGION_CLASS["full_ITS"]="./data/full_ITS/eukaryome_ITS_nocomplex.classification"
REGION_CLASS["ITS1"]="./data/ITS1/eukaryome_ITS1_nocomplex.classification"
REGION_CLASS["ITS2"]="./data/ITS2/eukaryome_ITS2_nocomplex.classification"

declare -A REGION_OUT
REGION_OUT["full_ITS"]="./data/full_ITS/eukaryome_cutoffs.txt"
REGION_OUT["ITS1"]="./data/ITS1/eukaryome_cutoffs.txt"
REGION_OUT["ITS2"]="./data/ITS2/eukaryome_cutoffs.txt"

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

echo "Activating conda environment..."
source ~/.bashrc
conda activate dyna_clust_predict

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if [[ ! -f "$CONSOLIDATE" ]]; then
    echo "ERROR: R script not found: $CONSOLIDATE" >&2
    exit 1
fi

# =============================================================================
# CONSOLIDATE EACH REGION
# =============================================================================

declare -a FAILED_REGIONS=()

for region in "${REGION_LABELS[@]}"; do

    cutoffs_in="${REGION_CUTOFFS_IN[$region]}"
    class_in="${REGION_CLASS[$region]}"
    out="${REGION_OUT[$region]}"

    echo ""
    echo "=== REGION: $region ==="
    echo "$(date)"

    if [[ ! -f "$cutoffs_in" ]]; then
        echo "WARNING: Cutoffs file not found, skipping region '$region': $cutoffs_in" >&2
        continue
    fi
    if [[ ! -f "$class_in" ]]; then
        echo "WARNING: Classification not found, skipping region '$region': $class_in" >&2
        continue
    fi

    Rscript "$CONSOLIDATE" \
        --cutoffs_in        "$cutoffs_in" \
        --classification_in "$class_in" \
        --output            "$out"

    rc=$?
    if [[ $rc -ne 0 ]]; then
        echo "ERROR: consolidate_cutoffs.R failed for region '$region' (exit code $rc)." >&2
        FAILED_REGIONS+=("$region")
    else
        echo "Finished region '$region' at: $(date)"
    fi

done

echo ""
if [[ ${#FAILED_REGIONS[@]} -gt 0 ]]; then
    echo "=== COMPLETED WITH ${#FAILED_REGIONS[@]} FAILURE(S): ${FAILED_REGIONS[*]} ==="
    exit 1
else
    echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="
fi
echo "$(date)"

conda deactivate
