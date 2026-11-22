#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=3950
#SBATCH --time=1-00:00:00
#SBATCH --partition=day
#SBATCH --output=logs/%x.%j.out

# Script name:  04_remove_complexes.sh
# Description:  For each ITS region, keep complete sequences (its_complete,
#               from 03) and remove species complexes: species the marker
#               cannot tell apart (R/remove_complexes.R).
#
#               Writes eukaryome_<region>_nocomplex.{fasta,classification},
#               which 05 to 08 read.
#
# Note:         This script must be run from the project root directory,
#               after 03 and before 05.

# =============================================================================
# PARAMETER SETUP
# =============================================================================

readonly REMOVE_COMPLEXES="./R/remove_complexes.R"
readonly N_CPUS="${SLURM_CPUS_PER_TASK:-$(nproc)}"

# Similarity at which two species are merged into a complex.
readonly THRESHOLD=1.0
# vsearch identity definition. 1 = terminal gaps count, so only identical
# sequences reach THRESHOLD=1.0. Same as 06, 07a and 07c.
readonly IDDEF=1
# Minimum overlap (bp) for a pair to count. Set to the minimum sequence length
# (02: 250 bp; 03: 50 bp). Same values as 06, 07a and 07c.
declare -A MIN_OVERLAP
MIN_OVERLAP["full_ITS"]=250
MIN_OVERLAP["ITS1"]=50
MIN_OVERLAP["ITS2"]=50
# A sequence matching this many species or more is left out of complex
# detection (likely mislabelled or contaminated).
readonly MAX_HUB_SPECIES=3
# Complex detection skips genera with fewer species than this.
readonly MIN_SPECIES_PER_PARENT=2
# Drop incomplete sequences (its_complete == FALSE) from the output.
readonly REQUIRE_COMPLETE="yes"
# Complex detection randomly downsamples larger genera to this many sequences.
readonly MAX_SEQ_NO=20000
readonly RUN_PARALLEL="yes"

REGION_LABELS=("full_ITS" "ITS1" "ITS2")

declare -A REGION_FASTA
REGION_FASTA["full_ITS"]="./data/full_ITS/eukaryome_ITS.fasta"
REGION_FASTA["ITS1"]="./data/ITS1/eukaryome_ITS1.fasta"
REGION_FASTA["ITS2"]="./data/ITS2/eukaryome_ITS2.fasta"

declare -A REGION_CLASS
REGION_CLASS["full_ITS"]="./data/full_ITS/eukaryome_ITS.classification"
REGION_CLASS["ITS1"]="./data/ITS1/eukaryome_ITS1.classification"
REGION_CLASS["ITS2"]="./data/ITS2/eukaryome_ITS2.classification"

declare -A REGION_FASTA_OUT
REGION_FASTA_OUT["full_ITS"]="./data/full_ITS/eukaryome_ITS_nocomplex.fasta"
REGION_FASTA_OUT["ITS1"]="./data/ITS1/eukaryome_ITS1_nocomplex.fasta"
REGION_FASTA_OUT["ITS2"]="./data/ITS2/eukaryome_ITS2_nocomplex.fasta"

declare -A REGION_CLASS_OUT
REGION_CLASS_OUT["full_ITS"]="./data/full_ITS/eukaryome_ITS_nocomplex.classification"
REGION_CLASS_OUT["ITS1"]="./data/ITS1/eukaryome_ITS1_nocomplex.classification"
REGION_CLASS_OUT["ITS2"]="./data/ITS2/eukaryome_ITS2_nocomplex.classification"

declare -A REGION_MANIFEST
REGION_MANIFEST["full_ITS"]="./data/full_ITS/complex_manifest.txt"
REGION_MANIFEST["ITS1"]="./data/ITS1/complex_manifest.txt"
REGION_MANIFEST["ITS2"]="./data/ITS2/complex_manifest.txt"

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

echo "Activating conda environment..."
source ~/.bashrc
conda activate dyna_clust_predict

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if [[ ! -f "$REMOVE_COMPLEXES" ]]; then
    echo "ERROR: R script not found: $REMOVE_COMPLEXES" >&2
    exit 1
fi

declare -a FAILED_REGIONS=()

# =============================================================================
# PROCESS EACH REGION
# =============================================================================

for region in "${REGION_LABELS[@]}"; do

    fasta_in="${REGION_FASTA[$region]}"
    class_in="${REGION_CLASS[$region]}"

    echo ""
    echo "=== REGION: $region ==="
    echo "$(date)"

    if [[ ! -f "$fasta_in" ]]; then
        echo "WARNING: FASTA not found, skipping region '$region': $fasta_in" >&2
        continue
    fi
    if [[ ! -f "$class_in" ]]; then
        echo "WARNING: Classification not found, skipping region '$region': $class_in" >&2
        continue
    fi

    Rscript "$REMOVE_COMPLEXES" \
        --fasta_in              "$fasta_in" \
        --classification_in     "$class_in" \
        --fasta_out              "${REGION_FASTA_OUT[$region]}" \
        --classification_out     "${REGION_CLASS_OUT[$region]}" \
        --manifest_out            "${REGION_MANIFEST[$region]}" \
        --threshold              "$THRESHOLD" \
        --iddef                  "$IDDEF" \
        --min_overlap            "${MIN_OVERLAP[$region]}" \
        --max_hub_species        "$MAX_HUB_SPECIES" \
        --min_species_per_parent "$MIN_SPECIES_PER_PARENT" \
        --require_complete       "$REQUIRE_COMPLETE" \
        --max_seq_no             "$MAX_SEQ_NO" \
        --n_cpus                 "$N_CPUS" \
        --run_parallel            "$RUN_PARALLEL" \
        --tmp_dir                "./tmp/remove_complexes_${region}"

    rc=$?
    if [[ $rc -ne 0 ]]; then
        echo "ERROR: remove_complexes.R failed for region '$region' (exit code $rc)." >&2
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
echo ""
echo "Output written to data/full_ITS, data/ITS1, data/ITS2"
echo "  Complex-and-completeness-filtered FASTA : eukaryome_<region>_nocomplex.fasta"
echo "  Matching classification                 : eukaryome_<region>_nocomplex.classification"
echo "  Manifest                                : complex_manifest.txt"

conda deactivate
