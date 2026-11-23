#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80
#SBATCH --time=7-00:00:00
#SBATCH --partition=week
#SBATCH --output=logs/%x.%j.out

# Script name:  07a_predict_cutoffs.sh
# Description:  Predict similarity cutoffs (R/predict.R) for full_ITS, ITS1 and
#               ITS2 in turn. For each region, runs every target x parent rank
#               combination (local), then each target rank with no parent
#               (global). Each dataset's FASTA and classification are built in
#               tmp/ from 05's ID files and deleted after use.
# Note:         This script must be run from the project root directory.

# =============================================================================
# PARAMETER SETUP
# =============================================================================

readonly PREDICT="./R/predict.R"
readonly TMP_DIR="./tmp"
readonly PREFIX="eukaryome"
# Threshold sweep: START_THRESH (per rank, below) to END_THRESH in steps of STEP
readonly STEP=0.001
readonly END_THRESH=1
readonly N_CPUS="${SLURM_CPUS_PER_TASK:-$(nproc)}"
readonly RUN_PARALLEL="yes"
# Minimum taxa with 2 or more sequences for a cutoff to be reported. A
# singleton scores a free Dice = 1.0 at threshold 1.0, so a dataset with few
# multi-sequence taxa has its optimum pinned near 1.0.
readonly MIN_MULTISEQ_GROUPS=10
# Round the reported cutoff to the nearest multiple of this (0.979 -> 0.98).
# The confidence is the F-measure at the rounded cutoff.
readonly CUTOFF_ROUND_TO=0.005
# vsearch identity definition. 1 = terminal gaps count, so a sequence that
# covers only part of another scores below 100%. Same as 04 and 06.
readonly IDDEF=1
# Recompute F-measures instead of reusing those cached in <prefix>.predicted;
# the cache is not refreshed when inputs or settings change.
readonly REDO="yes"

# Start thresholds per target rank
declare -A START_THRESH
START_THRESH["species"]="0.8"
START_THRESH["genus"]="0.7"
START_THRESH["family"]="0.6"
START_THRESH["order"]="0.5"
START_THRESH["class"]="0.5"
START_THRESH["phylum"]="0.5"

# Parent rank combinations per target rank
declare -A PARENT_RANKS
PARENT_RANKS["species"]="genus family order class phylum kingdom"
PARENT_RANKS["genus"]="family order class phylum kingdom"
PARENT_RANKS["family"]="order class phylum kingdom"
PARENT_RANKS["order"]="class phylum kingdom"
PARENT_RANKS["class"]="phylum kingdom"
PARENT_RANKS["phylum"]="kingdom"

# Filename abbreviations matching subset.R output names
declare -A RANK_ABBR
RANK_ABBR["kingdom"]="kng"
RANK_ABBR["phylum"]="phy"
RANK_ABBR["class"]="cls"
RANK_ABBR["order"]="ord"
RANK_ABBR["family"]="fam"
RANK_ABBR["genus"]="gen"

TARGET_RANKS=("species" "genus" "family" "order" "class" "phylum")
REGIONS=("full_ITS" "ITS1" "ITS2")

# Input FASTA per region (from 04)
declare -A REGION_FASTA
REGION_FASTA["full_ITS"]="./data/full_ITS/eukaryome_ITS_nocomplex.fasta"
REGION_FASTA["ITS1"]="./data/ITS1/eukaryome_ITS1_nocomplex.fasta"
REGION_FASTA["ITS2"]="./data/ITS2/eukaryome_ITS2_nocomplex.fasta"

# Input classification per region (from 04)
declare -A REGION_CLASS
REGION_CLASS["full_ITS"]="./data/full_ITS/eukaryome_ITS_nocomplex.classification"
REGION_CLASS["ITS1"]="./data/ITS1/eukaryome_ITS1_nocomplex.classification"
REGION_CLASS["ITS2"]="./data/ITS2/eukaryome_ITS2_nocomplex.classification"

# report_floor_hits <raw_cutoffs_file>
# Lists datasets whose cutoff is within 0.005 of START_THRESH; their optimum
# may lie below the sweep. Reports only: lowering START_THRESH can also move
# cutoffs where F is tied across a range, so it is a judgement call.
report_floor_hits() {
    local raw="$1"
    [[ -f "$raw" ]] || return 0
    echo ""
    echo "--- Cut-offs at the bottom of the sweep (<= start threshold + 0.005): $raw ---"
    awk -F'\t' -v sp="${START_THRESH[species]}" -v ge="${START_THRESH[genus]}" \
               -v fa="${START_THRESH[family]}" -v od="${START_THRESH[order]}" \
               -v cl="${START_THRESH[class]}" -v ph="${START_THRESH[phylum]}" '
        BEGIN { st["species"] = sp; st["genus"] = ge; st["family"] = fa
                st["order"] = od; st["class"] = cl; st["phylum"] = ph }
        NR == 1 { next }
        ($1 in st) && ($4 + 0 <= st[$1] + 0.005 + 1e-9) {
            n[$1]++; total++
            if (shown[$1]++ < 5) printf "  %-8s %-8s %-30s cut-off %s  F %s\n", $1, $2, $3, $4, $5
        }
        END {
            if (!total) { print "  none"; exit }
            for (r in n) printf "  %s: %d dataset(s)\n", r, n[r]
            print "  If many, consider a lower START_THRESH for that rank and re-run it."
        }' "$raw"
}

# =============================================================================
# DIRECTORY SETUP
# =============================================================================

mkdir -p "$TMP_DIR"

# =============================================================================
# HELPER FUNCTION
# =============================================================================

# subset_fasta_and_classification <id_file> <global_fasta> <global_class>
#                                 <out_fasta> <out_class>
# Subsets the region's FASTA and classification to the IDs in id_file.
subset_fasta_and_classification() {
    local id_file="$1"
    local global_fasta="$2"
    local global_class="$3"
    local out_fasta="$4"
    local out_class="$5"

    if [[ ! -f "$id_file" ]]; then
        echo "  WARNING: ID file not found, skipping: $id_file" >&2
        return 1
    fi

    # Keep sequences whose ID is in id_file
    seqkit grep -f "$id_file" "$global_fasta" > "$out_fasta"

    local n_seqs
    n_seqs=$(grep -c "^>" "$out_fasta")
    echo "  Wrote ${n_seqs} sequences to ${out_fasta}"

    # Keep the header row and rows whose ID is in id_file
    awk 'NR == FNR { ids[$1] = 1; next }
         FNR == 1  { print; next }
         $1 in ids { print }' "$id_file" "$global_class" > "$out_class"

    local n_class
    n_class=$(( $(wc -l < "$out_class") - 1 ))
    echo "  Wrote ${n_class} classification rows to ${out_class}"
}

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

echo "Activating conda environment..."
source ~/.bashrc
conda activate dyna_clust_predict

# Track failures for the end-of-job summary
declare -a FAILED_STEPS=()

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if [[ ! -f "$PREDICT" ]]; then
    echo "ERROR: R script not found: $PREDICT" >&2
    exit 1
fi

# =============================================================================
# MAIN LOOP: iterate over all three regions
# =============================================================================

for region in "${REGIONS[@]}"; do

    global_fasta="${REGION_FASTA[$region]}"
    global_class="${REGION_CLASS[$region]}"

    echo ""
    echo "=== REGION: $region ==="
    echo "$(date)"

    if [[ ! -f "$global_fasta" ]]; then
        echo "WARNING: FASTA not found, skipping region '$region': $global_fasta" >&2
        continue
    fi

    if [[ ! -f "$global_class" ]]; then
        echo "WARNING: Classification not found, skipping region '$region': $global_class" >&2
        continue
    fi

    # -------------------------------------------------------------------------
    # LOCAL PREDICTIONS (per target x parent combination)
    # -------------------------------------------------------------------------

    echo ""
    echo "--- LOCAL PREDICTIONS: $region ---"
    echo "$(date)"

    for target in "${TARGET_RANKS[@]}"; do

        st="${START_THRESH[$target]}"

        for parent in ${PARENT_RANKS[$target]}; do

            parent_abbr="${RANK_ABBR[$parent]}"
            id_file="./data/${region}/${target}_pred_id_${parent_abbr}.txt"
            tmp_fasta="$TMP_DIR/${target}_pred_${parent_abbr}.fasta"
            tmp_class="$TMP_DIR/${target}_pred_${parent_abbr}.classification"

            # No ID file: no dataset passed 05's filters
            if [[ ! -f "$id_file" ]]; then
                echo "  Skipping ${target} within ${parent}: ID file not found ($id_file)"
                continue
            fi

            echo ""
            echo "  -- ${region} / ${target} within ${parent} --"
            echo "  Subsetting sequences at: $(date)"
            subset_fasta_and_classification \
                "$id_file" "$global_fasta" "$global_class" \
                "$tmp_fasta" "$tmp_class" || continue

            echo "  Predicting cutoffs at: $(date)"
            Rscript "$PREDICT" \
                --input          "$tmp_fasta" \
                --classification "$tmp_class" \
                --rank           "$target" \
                --higher_rank    "$parent" \
                --start_threshold "${st}" \
                --end_threshold  "$END_THRESH" \
                --step           "$STEP" \
                --min_multiseq_groups "$MIN_MULTISEQ_GROUPS" \
                --iddef          "$IDDEF" \
                --redo           "$REDO" \
                --cutoff_round_to "$CUTOFF_ROUND_TO" \
                --prefix         "$PREFIX" \
                --id_col         id \
                --run_parallel   "$RUN_PARALLEL" \
                --n_cpus         "$N_CPUS" \
                --tmp_dir        "$TMP_DIR" \
                --out            "./data/${region}"

            rc=$?
            if [[ $rc -ne 0 ]]; then
                echo "" >&2
                echo "  ============================================================" >&2
                echo "  FAILED: ${region}/${target} within ${parent}  (exit code $rc)" >&2
                echo "  Time of failure: $(date)" >&2
                echo "  Memory at failure (RSS peak of this job):" >&2
                echo "    $(grep VmPeak /proc/$$/status 2>/dev/null || echo 'N/A')" >&2
                echo "    $(grep VmRSS  /proc/$$/status 2>/dev/null || echo 'N/A')" >&2
                if command -v sstat &>/dev/null && [[ -n "${SLURM_JOB_ID:-}" ]]; then
                    echo "  SLURM memory usage:" >&2
                    sstat -j "${SLURM_JOB_ID}.batch" --format=JobID,MaxRSS,MaxVMSize 2>/dev/null || true
                fi
                echo "  If the R output above mentions 'FutureInterruptError' or" >&2
                echo "  'parallel job did not deliver a result', the likely cause" >&2
                echo "  is the OS OOM killer terminating worker processes." >&2
                echo "  ============================================================" >&2
                echo "" >&2
                FAILED_STEPS+=("${region}/${target} within ${parent}")
            else
                echo "  Finished at: $(date)"
            fi

            rm -f "$tmp_fasta" "$tmp_class"

        done
    done

    # -------------------------------------------------------------------------
    # GLOBAL PREDICTIONS (no higher rank, per target rank)
    # -------------------------------------------------------------------------

    echo ""
    echo "--- GLOBAL PREDICTIONS: $region ---"
    echo "$(date)"

    for target in "${TARGET_RANKS[@]}"; do

        st="${START_THRESH[$target]}"
        id_file="./data/${region}/${target}_pred_id_global.txt"
        tmp_fasta="$TMP_DIR/${target}_unique.fasta"
        tmp_class="$TMP_DIR/${target}_unique.classification"

        if [[ ! -f "$id_file" ]]; then
            echo "  Skipping global prediction for ${target}: ID file not found ($id_file)"
            continue
        fi

        echo ""
        echo "  -- ${region} / global: ${target} --"
        echo "  Subsetting sequences at: $(date)"
        subset_fasta_and_classification \
            "$id_file" "$global_fasta" "$global_class" \
            "$tmp_fasta" "$tmp_class" || continue

        echo "  Predicting global cutoffs at: $(date)"
        Rscript "$PREDICT" \
            --input          "$tmp_fasta" \
            --classification "$tmp_class" \
            --rank           "$target" \
            --start_threshold "${st}" \
            --end_threshold  "$END_THRESH" \
            --step           "$STEP" \
            --min_multiseq_groups "$MIN_MULTISEQ_GROUPS" \
            --iddef          "$IDDEF" \
            --redo           "$REDO" \
            --cutoff_round_to "$CUTOFF_ROUND_TO" \
            --prefix         "$PREFIX" \
            --id_col         id \
            --run_parallel   "$RUN_PARALLEL" \
            --n_cpus         "$N_CPUS" \
            --tmp_dir        "$TMP_DIR" \
            --out            "./data/${region}"

        rc=$?
        if [[ $rc -ne 0 ]]; then
            echo "" >&2
            echo "  ============================================================" >&2
            echo "  FAILED: ${region}/global ${target}  (exit code $rc)" >&2
            echo "  Time of failure: $(date)" >&2
            echo "  ============================================================" >&2
            echo "" >&2
            FAILED_STEPS+=("${region}/global ${target}")
        else
            echo "  Finished at: $(date)"
        fi

        rm -f "$tmp_fasta" "$tmp_class"

    done

    report_floor_hits "./data/${region}/${PREFIX}_raw_cutoffs.txt"

    echo ""
    echo "=== REGION $region COMPLETE at $(date) ==="

done

# =============================================================================
# CLEANUP
# =============================================================================

echo ""
echo "=== CLEANUP ==="
echo "Removing tmp directory..."
rm -rf "$TMP_DIR"

echo ""
if [[ ${#FAILED_STEPS[@]} -gt 0 ]]; then
    echo "=== PIPELINE COMPLETED WITH ${#FAILED_STEPS[@]} FAILURE(S) ==="
    echo "$(date)"
    echo ""
    echo "Failed steps:"
    for step in "${FAILED_STEPS[@]}"; do
        echo "  - ${step}"
    done
    echo ""
    echo "Likely cause: OS OOM (out-of-memory) killer terminated worker processes."
    echo "Possible fixes:"
    echo "  1. Request more memory in SLURM (#SBATCH --mem=...)."
    echo "  2. Lower MAX_SEQUENCES in 05_prepare_subsets.sh and re-run it."
    echo "  3. Reduce --n_cpus to run fewer parallel workers."
else
    echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="
    echo "$(date)"
fi
echo ""
echo "Cutoff files written to data/full_ITS, data/ITS1, data/ITS2"
echo "  Local cutoffs  : ${PREFIX}_raw_cutoffs.json (per region)"
echo "  Global cutoffs : ${PREFIX}_raw_cutoffs.json (appended per region)"

conda deactivate
