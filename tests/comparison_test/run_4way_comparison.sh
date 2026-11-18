#!/bin/bash
# run_4way_comparison.sh
# Compare 4 sequence selection approaches on prediction results:
#   1. random:          Pure random sample
#   2. random_capped:   Random + dominant group capped at ≤50%
#   3. structured:      Hierarchical budget allocation
#   4. structured_min5: Hierarchical + each subgroup ≥5 seqs
set -euo pipefail

PROJECT_DIR="/Users/lukelikesdirt/Desktop/project/dyna-clust-predict"
DATA_DIR="${PROJECT_DIR}/tests/comparison_test/data"
FASTA="${DATA_DIR}/test_selection.fasta"
CLASS="${DATA_DIR}/test_selection.classification"
RANKS=("genus" "family" "order" "class" "phylum")
APPROACHES=("random" "random_capped" "structured" "structured_min5")

cd "$PROJECT_DIR"

# Create result directories
for APP in "${APPROACHES[@]}"; do
    mkdir -p "${PROJECT_DIR}/tests/comparison_test/results_${APP}"
done

for RANK in "${RANKS[@]}"; do
    echo "========================================"
    echo "RANK: $RANK"
    echo "========================================"

    for APP in "${APPROACHES[@]}"; do
        IDS="${DATA_DIR}/sel4_${APP}_${RANK}_ids.txt"
        APP_FASTA="${DATA_DIR}/sel4_${APP}_${RANK}.fasta"
        APP_SIM="${DATA_DIR}/sel4_${APP}_${RANK}.sim"
        RESULT_DIR="${PROJECT_DIR}/tests/comparison_test/results_${APP}"

        if [[ ! -f "$IDS" ]]; then
            echo "  [$APP] SKIP: no ID file"
            continue
        fi

        N_IDS=$(wc -l < "$IDS" | tr -d ' ')
        if [[ "$N_IDS" -lt 30 ]]; then
            echo "  [$APP] SKIP: only $N_IDS IDs (need ≥30)"
            continue
        fi

        # Extract FASTA subset
        if [[ ! -f "$APP_FASTA" ]]; then
            echo "  [$APP] Extracting FASTA subset..."
            seqkit grep -f "$IDS" "$FASTA" > "$APP_FASTA"
        fi
        N_SEQ=$(grep -c '^>' "$APP_FASTA")
        echo "  [$APP] $N_SEQ seqs"

        # Compute sim matrix
        if [[ ! -f "$APP_SIM" ]]; then
            echo "  [$APP] Computing sim matrix..."
            Rscript R/compute_sim.R -i "$APP_FASTA" -o "$DATA_DIR" --n_cpus 4 2>&1 | tail -1
        fi

        # Run prediction
        echo "  [$APP] Predicting..."
        Rscript R/predict.R \
            -i "$APP_FASTA" \
            -c "$CLASS" \
            --rank "$RANK" \
            -o "$RESULT_DIR" \
            --prefix "sel_${RANK}" \
            --sim "$APP_SIM" \
            --start_threshold 0 --end_threshold 1 --step 0.01 \
            --min_group_no 10 --min_seq_no 30 --max_proportion 0.5 \
            --run_parallel no --n_cpus 4 \
            2>&1 | grep -E "(cutoff=|Rank |Given )" || true

        echo ""
    done
done

# ── Summary table ─────────────────────────────────────────────────────────────
echo ""
echo "═══════════════════════════════════════════════════════════════════════"
echo "                        4-WAY COMPARISON SUMMARY"
echo "═══════════════════════════════════════════════════════════════════════"
printf "%-10s  %-18s  %8s  %10s  %6s  %8s  %8s\n" \
    "Rank" "Approach" "Cutoff" "Confidence" "Seqs" "Groups" "MaxProp"
echo "───────────────────────────────────────────────────────────────────────"

for RANK in "${RANKS[@]}"; do
    for APP in "${APPROACHES[@]}"; do
        RESULT_FILE="${PROJECT_DIR}/tests/comparison_test/results_${APP}/sel_${RANK}.cutoffs.json.txt"
        if [[ -f "$RESULT_FILE" ]]; then
            tail -n +2 "$RESULT_FILE" | while IFS=$'\t' read -r _ _ _ cutoff conf seqs groups maxp _; do
                printf "%-10s  %-18s  %8s  %10s  %6s  %8s  %8s\n" \
                    "$RANK" "$APP" "$cutoff" "$conf" "$seqs" "$groups" "$maxp"
            done
        else
            printf "%-10s  %-18s  %8s  %10s  %6s  %8s  %8s\n" \
                "$RANK" "$APP" "-" "-" "-" "-" "-"
        fi
    done
    echo "───────────────────────────────────────────────────────────────────────"
done
