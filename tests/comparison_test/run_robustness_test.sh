#!/bin/bash
# run_robustness_test.sh
# Run the 4-way selection comparison across multiple independent test datasets
set -euo pipefail

PROJECT_DIR="/Users/lukelikesdirt/Desktop/project/dyna-clust-predict"
ROBUST_DIR="${PROJECT_DIR}/tests/comparison_test/robustness"
RANKS=("genus" "family" "order" "class" "phylum")
APPROACHES=("random" "random_capped" "structured" "structured_min5")
MAX_SEQS=500

cd "$PROJECT_DIR"

N_DS=$(ls -d "${ROBUST_DIR}"/dataset_* 2>/dev/null | wc -l | tr -d ' ')
echo "Running 4-way comparison across $N_DS datasets"
echo ""

for DS_DIR in "${ROBUST_DIR}"/dataset_*; do
    DS_NAME=$(basename "$DS_DIR")
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║  $DS_NAME"
    echo "╚══════════════════════════════════════════════════════════════╝"

    FASTA="${DS_DIR}/test.fasta"
    CLASS="${DS_DIR}/test.classification"

    # Step 1: Generate ID lists for all 4 approaches
    echo "  Generating selection ID lists..."
    Rscript tests/comparison_test/selection_test_4way.R \
        "$FASTA" "$CLASS" "$DS_DIR" "$MAX_SEQS" 2>&1 | sed 's/^/    /'

    # Step 2: For each rank × approach, extract FASTA, compute sim, predict
    for RANK in "${RANKS[@]}"; do
        echo "  ── $RANK ──"
        for APP in "${APPROACHES[@]}"; do
            IDS="${DS_DIR}/sel4_${APP}_${RANK}_ids.txt"
            APP_FASTA="${DS_DIR}/sel4_${APP}_${RANK}.fasta"
            APP_SIM="${DS_DIR}/sel4_${APP}_${RANK}.sim"
            RES_DIR="${DS_DIR}/results_${APP}"
            mkdir -p "$RES_DIR"

            if [[ ! -f "$IDS" ]]; then
                echo "    [$APP] SKIP: no ID file"
                continue
            fi

            N_IDS=$(wc -l < "$IDS" | tr -d ' ')
            if [[ "$N_IDS" -lt 30 ]]; then
                echo "    [$APP] SKIP: only $N_IDS IDs"
                continue
            fi

            # Extract FASTA
            if [[ ! -f "$APP_FASTA" ]]; then
                seqkit grep -f "$IDS" "$FASTA" > "$APP_FASTA" 2>/dev/null
            fi

            # Compute sim
            if [[ ! -f "$APP_SIM" ]]; then
                Rscript R/compute_sim.R -i "$APP_FASTA" -o "$DS_DIR" --n_cpus 4 2>/dev/null
            fi

            # Predict
            Rscript R/predict.R \
                -i "$APP_FASTA" \
                -c "$CLASS" \
                --rank "$RANK" \
                -o "$RES_DIR" \
                --prefix "sel_${RANK}" \
                --sim "$APP_SIM" \
                --start_threshold 0 --end_threshold 1 --step 0.01 \
                --min_group_no 10 --min_seq_no 30 --max_proportion 0.5 \
                --run_parallel no --n_cpus 4 \
                2>/dev/null

            # Extract result
            RES_FILE="${RES_DIR}/sel_${RANK}.cutoffs.json.txt"
            if [[ -f "$RES_FILE" ]]; then
                RESULT=$(tail -n +2 "$RES_FILE" | awk -F'\t' '{printf "cutoff=%s conf=%s seqs=%s groups=%s", $4, $5, $6, $7}')
                echo "    [$APP] $RESULT"
            else
                echo "    [$APP] NO RESULT"
            fi
        done
    done
    echo ""
done

# ── Cross-dataset summary ────────────────────────────────────────────────────
echo ""
echo "╔══════════════════════════════════════════════════════════════════════════╗"
echo "║                    ROBUSTNESS SUMMARY (ALL DATASETS)                   ║"
echo "╚══════════════════════════════════════════════════════════════════════════╝"
echo ""

for RANK in "${RANKS[@]}"; do
    echo "── $RANK ──"
    printf "  %-18s" "Approach"
    for DS_DIR in "${ROBUST_DIR}"/dataset_*; do
        printf "  %10s" "$(basename "$DS_DIR" | sed 's/dataset_/ds/')"
    done
    printf "  %10s  %10s\n" "Mean" "SD"
    echo "  $(printf '%.0s─' {1..90})"

    for APP in "${APPROACHES[@]}"; do
        printf "  %-18s" "$APP"
        CUTOFFS=()
        CONFS=()
        for DS_DIR in "${ROBUST_DIR}"/dataset_*; do
            RES_FILE="${DS_DIR}/results_${APP}/sel_${RANK}.cutoffs.json.txt"
            if [[ -f "$RES_FILE" ]]; then
                CUTOFF=$(tail -n +2 "$RES_FILE" | awk -F'\t' '{print $4}')
                CONF=$(tail -n +2 "$RES_FILE" | awk -F'\t' '{print $5}')
                printf "  %4s/%5s" "$CUTOFF" "$CONF"
                CUTOFFS+=("$CUTOFF")
                CONFS+=("$CONF")
            else
                printf "  %10s" "-"
            fi
        done

        # Compute mean/SD for cutoffs
        if [[ ${#CUTOFFS[@]} -gt 0 ]]; then
            MEAN_CUT=$(printf '%s\n' "${CUTOFFS[@]}" | awk '{s+=$1} END {printf "%.2f", s/NR}')
            SD_CUT=$(printf '%s\n' "${CUTOFFS[@]}" | awk -v m="$MEAN_CUT" '{d=$1-m; s+=d*d} END {printf "%.3f", sqrt(s/NR)}')
            MEAN_CONF=$(printf '%s\n' "${CONFS[@]}" | awk '{s+=$1} END {printf "%.3f", s/NR}')
            SD_CONF=$(printf '%s\n' "${CONFS[@]}" | awk -v m="$MEAN_CONF" '{d=$1-m; s+=d*d} END {printf "%.3f", sqrt(s/NR)}')
            printf "  %4s/%5s  %5s/%5s" "$MEAN_CUT" "$MEAN_CONF" "$SD_CUT" "$SD_CONF"
        fi
        echo ""
    done
    echo ""
done
