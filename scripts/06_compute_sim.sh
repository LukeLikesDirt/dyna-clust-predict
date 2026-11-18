#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80
#SBATCH --time=7-00:00:00
#SBATCH --partition=week
#SBATCH --output=logs/%x.%j.out

# Script name:  06_compute_sim.sh
# Description:  Optional. Pre-compute pairwise vsearch similarity (.sim) files
#               from each rank's unique-sequence IDs (05 STEP 1), for full_ITS,
#               ITS1 and ITS2 (R/compute_sim.R). Ranks with an existing .sim
#               file are skipped. Pass a .sim file to predict.R with --sim to
#               skip vsearch there (07a and 07c do not).
#
# Note:         This script must be run from the project root directory.

# =============================================================================
# PARAMETER SETUP
# =============================================================================

readonly SIM="./R/compute_sim.R"
readonly TMP_DIR="./tmp"
readonly MIN_SIM=0        # lowest similarity written to the .sim file
readonly IDDEF=2          # vsearch identity definition (terminal gaps excluded); same as 04, 07a and 07c
readonly N_CPUS="${SLURM_CPUS_PER_TASK:-$(nproc)}"

# Minimum overlap (bp) for a pair to be written; without it, sequences that
# share a few bases score 100% under IDDEF=2. Same values as 04, 07a and 07c.
declare -A MIN_OVERLAP
MIN_OVERLAP["full_ITS"]=250
MIN_OVERLAP["ITS1"]=50
MIN_OVERLAP["ITS2"]=50

TARGET_RANKS=("species" "genus" "family" "order" "class" "phylum")
REGIONS=("full_ITS" "ITS1" "ITS2")

# =============================================================================
# DIRECTORY SETUP
# =============================================================================

mkdir -p "$TMP_DIR"

# =============================================================================
# HELPER FUNCTION: resolve FASTA path for a given region
# =============================================================================

region_fasta() {
    local region="$1"
    case "$region" in
        full_ITS) echo "./data/full_ITS/eukaryome_ITS_nocomplex.fasta"  ;;
        ITS1)     echo "./data/ITS1/eukaryome_ITS1_nocomplex.fasta"     ;;
        ITS2)     echo "./data/ITS2/eukaryome_ITS2_nocomplex.fasta"     ;;
        *)
            echo "ERROR: Unknown region: $region" >&2
            return 1
            ;;
    esac
}

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

echo "Activating conda environment..."
source ~/.bashrc
conda activate dyna_clust_predict

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if [[ ! -f "$SIM" ]]; then
    echo "ERROR: R script not found: $SIM" >&2
    exit 1
fi

# =============================================================================
# COMPUTE SIMILARITY MATRICES
# =============================================================================

echo ""
echo "=== COMPUTING VSEARCH SIMILARITY MATRICES ==="
echo "$(date)"
echo "Threads : $N_CPUS"
echo "min_sim : $MIN_SIM"
echo "iddef   : $IDDEF"
echo "min_overlap (bp): full_ITS=${MIN_OVERLAP[full_ITS]} ITS1=${MIN_OVERLAP[ITS1]} ITS2=${MIN_OVERLAP[ITS2]}"
echo ""

for region in "${REGIONS[@]}"; do

    global_fasta=$(region_fasta "$region")

    echo "=== REGION: $region ==="

    if [[ ! -f "$global_fasta" ]]; then
        echo "  WARNING: FASTA not found, skipping region '$region': $global_fasta" >&2
        echo ""
        continue
    fi

    for target in "${TARGET_RANKS[@]}"; do

        id_file="./data/${region}/${target}_unique_id.txt"
        sim_out="./data/${region}/${target}_unique.sim"
        tmp_fasta="$TMP_DIR/${region}_${target}_unique.fasta"

        echo "--- ${region} / ${target} ---"

        # Skip if sim file already exists
        if [[ -f "$sim_out" ]]; then
            echo "  Sim file already exists, skipping: $sim_out"
            continue
        fi

        # Skip if ID file is missing
        if [[ ! -f "$id_file" ]]; then
            echo "  Skipping: ID file not found ($id_file)"
            continue
        fi

        echo "  Subsetting sequences at: $(date)"
        seqkit grep -f "$id_file" "$global_fasta" > "$tmp_fasta"

        n_seqs=$(grep -c "^>" "$tmp_fasta")
        echo "  Sequences in subset: $n_seqs"

        echo "  Computing similarity at: $(date)"
        Rscript "$SIM" \
            --input  "$tmp_fasta" \
            --out    "./data/${region}" \
            --min_sim "$MIN_SIM" \
            --iddef   "$IDDEF" \
            --min_overlap "${MIN_OVERLAP[$region]}" \
            --n_cpus  "$N_CPUS"

        if [[ $? -ne 0 ]]; then
            echo "  WARNING: compute_sim.R failed for ${region}/${target}" >&2
        else
            # compute_sim.R names its output after the tmp FASTA; rename it.
            local_sim="./data/${region}/${region}_${target}_unique.sim"
            if [[ -f "$local_sim" && "$local_sim" != "$sim_out" ]]; then
                mv "$local_sim" "$sim_out"
            fi
            echo "  Finished at: $(date)"
            echo "  Sim file: $sim_out"
        fi

        rm -f "$tmp_fasta"
        echo ""

    done

    echo ""

done

# =============================================================================
# CLEANUP
# =============================================================================

echo "=== CLEANUP ==="
echo "Removing tmp directory..."
rm -rf "$TMP_DIR"

echo ""
echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="
echo "$(date)"
echo ""
echo "Sim files written per region:"
for region in "${REGIONS[@]}"; do
    for target in "${TARGET_RANKS[@]}"; do
        sim_out="./data/${region}/${target}_unique.sim"
        [[ -f "$sim_out" ]] && echo "  $sim_out"
    done
done

conda deactivate
