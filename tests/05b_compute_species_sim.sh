#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80
#SBATCH --time=60-00:00:00
#SBATCH --partition=long
#SBATCH --output=logs/%x.%j.out

# Script name:  05b_compute_species_sim.sh
# Description:  Pre-compute a pairwise vsearch similarity matrix for all
#               sequences identified at species rank, for each ITS region.
#               Species-identified IDs are extracted directly from the
#               classification file (matching the is_identified logic in utils.R)
#               and used to subset the full FASTA before running compute_sim.R.
#
#               Outputs (one per region):
#                 ./data/full_ITS/species_ITS.sim
#                 ./data/ITS1/species_ITS1.sim
#                 ./data/ITS2/species_ITS2.sim
#
#               Skips any region whose .sim file already exists.
# Note:         This script must be run from the project root directory.

# =============================================================================
# PARAMETER SETUP
# =============================================================================

readonly SIM="./R/compute_sim.R"
readonly TMP_DIR="./tmp"
readonly MIN_SIM=0
readonly N_CPUS="${SLURM_CPUS_PER_TASK:-$(nproc)}"

# Parallel arrays: label, input FASTA, classification file, output sim name
REGION_LABELS=("full_ITS"                              "ITS1"                          "ITS2")
REGION_FASTAS=("./data/full_ITS/eukaryome_ITS.fasta"   "./data/ITS1/eukaryome_ITS1.fasta"  "./data/ITS2/eukaryome_ITS2.fasta")
REGION_CLASSIFICATIONS=("./data/full_ITS/eukaryome_ITS.classification" \
                        "./data/ITS1/eukaryome_ITS1.classification" \
                        "./data/ITS2/eukaryome_ITS2.classification")
REGION_SIMNAMES=("species_ITS" "species_ITS1" "species_ITS2")

# =============================================================================
# DIRECTORY SETUP
# =============================================================================

mkdir -p "$TMP_DIR"

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
# HELPER: extract species-identified IDs from a classification file
# Mirrors the is_identified() logic from R/utils.R (case-insensitive).
# Expects a tab-delimited file with a header row containing an "id" column
# and a "species" column.
# =============================================================================

extract_species_ids() {
    local classification_file="$1"
    awk -F'\t' '
    NR == 1 {
        for (i = 1; i <= NF; i++) col[$i] = i
        id_col  = col["id"]
        sp_col  = col["species"]
        next
    }
    {
        sp = $sp_col
        if (sp == "" || sp == "NA") next
        sl = tolower(sp)
        if (sl ~ /^unidentified$/ || sl ~ /^unclassified$/) next
        if (sl ~ /^uncultured/)                              next
        if (sl ~ /incertae[ _.-]sedis/)                     next
        if (sp ~ /[ _.-]sp\./)                              next
        if (sp ~ /[Uu]nispike[123]/)                        next
        if (sl ~ /archaea|bacteria/)                        next
        if (sl ~ /mitochondrion|nucleomorph|plastid/)       next
        print $id_col
    }' "$classification_file"
}

# =============================================================================
# COMPUTE SIMILARITY MATRICES
# =============================================================================

echo ""
echo "=== COMPUTING VSEARCH SPECIES SIMILARITY MATRICES ==="
echo "$(date)"
echo "Threads : $N_CPUS"
echo "min_sim : $MIN_SIM"
echo ""

for i in "${!REGION_LABELS[@]}"; do

    label="${REGION_LABELS[$i]}"
    fasta="${REGION_FASTAS[$i]}"
    classification="${REGION_CLASSIFICATIONS[$i]}"
    sim_name="${REGION_SIMNAMES[$i]}"
    out_dir="./data/${label}"
    sim_out="${out_dir}/${sim_name}.sim"
    tmp_fasta="${TMP_DIR}/${sim_name}.fasta"

    echo "=== REGION: $label ==="

    # Skip if sim file already exists
    if [[ -f "$sim_out" ]]; then
        echo "  Sim file already exists, skipping: $sim_out"
        echo ""
        continue
    fi

    if [[ ! -f "$fasta" ]]; then
        echo "  WARNING: FASTA not found, skipping: $fasta" >&2
        echo ""
        continue
    fi

    if [[ ! -f "$classification" ]]; then
        echo "  WARNING: Classification not found, skipping: $classification" >&2
        echo ""
        continue
    fi

    echo "  Extracting species-identified IDs at: $(date)"
    id_file="${TMP_DIR}/${sim_name}_ids.txt"
    extract_species_ids "$classification" > "$id_file"

    n_ids=$(wc -l < "$id_file")
    echo "  Species-identified sequences: $n_ids"

    if [[ "$n_ids" -eq 0 ]]; then
        echo "  WARNING: No species-identified sequences found, skipping." >&2
        rm -f "$id_file"
        echo ""
        continue
    fi

    echo "  Subsetting FASTA at: $(date)"
    seqkit grep -f "$id_file" "$fasta" > "$tmp_fasta"
    rm -f "$id_file"

    n_seqs=$(grep -c "^>" "$tmp_fasta" || true)
    echo "  Sequences in subset: $n_seqs"

    echo "  Computing similarity at: $(date)"
    Rscript "$SIM" \
        --input    "$tmp_fasta" \
        --out      "$out_dir" \
        --min_sim  "$MIN_SIM" \
        --n_cpus   "$N_CPUS"

    if [[ $? -ne 0 ]]; then
        echo "  WARNING: compute_sim.R failed for $label" >&2
    else
        echo "  Finished at: $(date)"
        echo "  Sim file: $sim_out"
    fi

    rm -f "$tmp_fasta"
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
echo "Sim files written:"
for i in "${!REGION_LABELS[@]}"; do
    sim_out="./data/${REGION_LABELS[$i]}/${REGION_SIMNAMES[$i]}.sim"
    [[ -f "$sim_out" ]] && echo "  $sim_out"
done

conda deactivate
