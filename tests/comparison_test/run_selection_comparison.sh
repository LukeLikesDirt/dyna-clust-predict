#!/bin/bash
# run_selection_comparison.sh
# Compare structured vs random selection on prediction results
set -euo pipefail

PROJECT_DIR="/Users/lukelikesdirt/Desktop/project/dyna-clust-predict"
DATA_DIR="${PROJECT_DIR}/tests/comparison_test/data"
RES_STRUCTURED="${PROJECT_DIR}/tests/comparison_test/results_structured"
RES_RANDOM="${PROJECT_DIR}/tests/comparison_test/results_random"
FASTA="${DATA_DIR}/test_selection.fasta"
CLASS="${DATA_DIR}/test_selection.classification"
RANKS=("genus" "family" "order" "class" "phylum")

cd "$PROJECT_DIR"
mkdir -p "$RES_STRUCTURED" "$RES_RANDOM"

for RANK in "${RANKS[@]}"; do
    STRUCTURED_IDS="${DATA_DIR}/selection_structured_${RANK}_ids.txt"
    RANDOM_IDS="${DATA_DIR}/selection_random_${RANK}_ids.txt"

    if [[ ! -f "$STRUCTURED_IDS" ]]; then
        echo "SKIP $RANK: no structured ID file"
        continue
    fi

    echo "========================================"
    echo "RANK: $RANK"
    echo "========================================"

    # Extract FASTA subsets
    STRUCT_FASTA="${DATA_DIR}/sel_struct_${RANK}.fasta"
    RANDOM_FASTA="${DATA_DIR}/sel_random_${RANK}.fasta"

    if [[ ! -f "$STRUCT_FASTA" ]]; then
        seqkit grep -f "$STRUCTURED_IDS" "$FASTA" > "$STRUCT_FASTA"
    fi
    if [[ ! -f "$RANDOM_FASTA" ]]; then
        seqkit grep -f "$RANDOM_IDS" "$FASTA" > "$RANDOM_FASTA"
    fi
    echo "  Structured FASTA: $(grep -c '^>' "$STRUCT_FASTA") seqs"
    echo "  Random FASTA:     $(grep -c '^>' "$RANDOM_FASTA") seqs"

    # Compute sim matrices
    STRUCT_SIM="${DATA_DIR}/sel_struct_${RANK}.sim"
    RANDOM_SIM="${DATA_DIR}/sel_random_${RANK}.sim"

    if [[ ! -f "$STRUCT_SIM" ]]; then
        echo "  Computing structured sim matrix..."
        Rscript R/compute_sim.R -i "$STRUCT_FASTA" -o "$DATA_DIR" --n_cpus 4 2>&1 | tail -1
        mv "${DATA_DIR}/sel_struct_${RANK}.sim" "$STRUCT_SIM" 2>/dev/null || true
    fi
    if [[ ! -f "$RANDOM_SIM" ]]; then
        echo "  Computing random sim matrix..."
        Rscript R/compute_sim.R -i "$RANDOM_FASTA" -o "$DATA_DIR" --n_cpus 4 2>&1 | tail -1
        mv "${DATA_DIR}/sel_random_${RANK}.sim" "$RANDOM_SIM" 2>/dev/null || true
    fi

    # Run predictions
    echo "  Predicting (structured)..."
    Rscript R/predict.R \
        -i "$STRUCT_FASTA" \
        -c "$CLASS" \
        --rank "$RANK" \
        -o "$RES_STRUCTURED" \
        --prefix "sel_${RANK}" \
        --sim "$STRUCT_SIM" \
        --start_threshold 0 --end_threshold 1 --step 0.01 \
        --min_group_no 10 --min_seq_no 30 --max_proportion 0.5 \
        --run_parallel no --n_cpus 4 \
        2>&1 | grep -E "(cutoff=|Rank )"

    echo "  Predicting (random)..."
    Rscript R/predict.R \
        -i "$RANDOM_FASTA" \
        -c "$CLASS" \
        --rank "$RANK" \
        -o "$RES_RANDOM" \
        --prefix "sel_${RANK}" \
        --sim "$RANDOM_SIM" \
        --start_threshold 0 --end_threshold 1 --step 0.01 \
        --min_group_no 10 --min_seq_no 30 --max_proportion 0.5 \
        --run_parallel no --n_cpus 4 \
        2>&1 | grep -E "(cutoff=|Rank )"

    echo ""
done

echo "========================================"
echo "SUMMARY"
echo "========================================"
echo ""
echo "Structured selection results:"
for f in "$RES_STRUCTURED"/sel_*.cutoffs.json.txt; do
    if [[ -f "$f" ]]; then
        echo "  $(basename "$f"):"
        tail -n +2 "$f" | awk -F'\t' '{printf "    cutoff=%s  confidence=%s  seqs=%s  groups=%s  maxprop=%s\n", $4, $5, $6, $7, $8}'
    fi
done

echo ""
echo "Random selection results:"
for f in "$RES_RANDOM"/sel_*.cutoffs.json.txt; do
    if [[ -f "$f" ]]; then
        echo "  $(basename "$f"):"
        tail -n +2 "$f" | awk -F'\t' '{printf "    cutoff=%s  confidence=%s  seqs=%s  groups=%s  maxprop=%s\n", $4, $5, $6, $7, $8}'
    fi
done
