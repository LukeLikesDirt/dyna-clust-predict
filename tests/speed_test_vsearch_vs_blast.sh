#!/usr/bin/env bash
# speed_test_vsearch_vs_blast.sh — Compare vsearch vs blastn-short approaches.
#
# Runs predict.R (vsearch) and predict_blast.R (blastn-short) on the same
# dataset(s) and reports wall-clock time, cutoff values, and F-measures.
#
# Usage:
#   bash tests/speed_test_vsearch_vs_blast.sh [small|large|both]
#
# Modes:
#   small  — 40 sequences  (tmp/test.fasta, global prediction)
#   large  — 620 sequences (tmp/bench/local_many.fasta, local prediction by genus)
#   both   — run both (default)
#
# Outputs: tmp/speed_test/ directory with timing and result files.

set -euo pipefail

cd "$(dirname "$0")/.."

MODE="${1:-both}"
OUTDIR="tmp/speed_test"
N_CPUS="${N_CPUS:-4}"
REPS="${REPS:-3}"

mkdir -p "$OUTDIR"

# ── Colour helpers ────────────────────────────────────────────────────────────
BOLD='\033[1m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
RESET='\033[0m'

log()  { echo -e "${CYAN}[speed_test]${RESET} $*"; }
bold() { echo -e "${BOLD}$*${RESET}"; }

# ── Run a single timed invocation ─────────────────────────────────────────────
# Usage: run_timed <label> <script> <args...>
# Captures wall time, writes output to $OUTDIR/<label>.log
run_timed() {
  local label="$1"; shift
  local script="$1"; shift
  local logfile="$OUTDIR/${label}.log"

  log "Running: $label"

  # GNU time on macOS is gtime (brew install gnu-time), fall back to /usr/bin/time
  local TIME_CMD
  if command -v gtime &>/dev/null; then
    TIME_CMD="gtime"
  else
    TIME_CMD="/usr/bin/time"
  fi

  # Capture wall-clock seconds
  local start_ns end_ns elapsed
  start_ns=$(python3 -c "import time; print(time.time())")
  Rscript "$script" "$@" > "$logfile" 2>&1 || true
  end_ns=$(python3 -c "import time; print(time.time())")
  elapsed=$(python3 -c "print(round($end_ns - $start_ns, 2))")

  echo "$elapsed" > "$OUTDIR/${label}.time"
  log "  → $label finished in ${elapsed}s"
}

# ── Extract cutoff from JSON ──────────────────────────────────────────────────
get_cutoff() {
  local json_file="$1"
  python3 -c "
import json, sys
with open('$json_file') as f:
    d = json.load(f)
for rank in d:
    for ds in d[rank]:
        v = d[rank][ds]
        co = v.get('cut-off', 'N/A')
        cf = v.get('confidence', 'N/A')
        print(f'{rank}/{ds}: cutoff={co}  F={cf}')
" 2>/dev/null || echo "(no results)"
}

# ── Small test: 40 sequences, global ─────────────────────────────────────────
run_small() {
  bold "\n════ Small dataset: 40 sequences (global) ════"

  local FASTA="tmp/test.fasta"
  local CLASS="tmp/test.classification"
  local RANK="species"
  local START_T="0.8"
  local END_T="1.0"

  local vsearch_times=()
  local blast_times=()

  for rep in $(seq 1 "$REPS"); do
    log "Repetition $rep/$REPS"

    # Clean outputs between runs
    rm -f "$OUTDIR"/small_vsearch_r${rep}.*
    rm -f "$OUTDIR"/small_blast_r${rep}.*

    run_timed "small_vsearch_r${rep}" R/predict.R \
      -i "$FASTA" -c "$CLASS" --rank "$RANK" \
      --start_threshold "$START_T" --end_threshold "$END_T" \
      --min_seq_no 2 --min_group_no 2 \
      -o "$OUTDIR" --prefix "small_vsearch_r${rep}" \
      --n_cpus "$N_CPUS" --run_parallel no --redo yes \
      --tmp_dir "$OUTDIR/tmp_vs"

    run_timed "small_blast_r${rep}" R/predict_blast.R \
      -i "$FASTA" -c "$CLASS" --rank "$RANK" \
      --start_threshold "$START_T" --end_threshold "$END_T" \
      --min_seq_no 2 --min_group_no 2 \
      -o "$OUTDIR" --prefix "small_blast_r${rep}" \
      --n_cpus "$N_CPUS" --run_parallel no --redo yes \
      --tmp_dir "$OUTDIR/tmp_bl"

    vsearch_times+=("$(cat "$OUTDIR/small_vsearch_r${rep}.time")")
    blast_times+=("$(cat "$OUTDIR/small_blast_r${rep}.time")")
  done

  bold "\n── Small dataset results ──"
  log "vsearch cutoffs:"
  get_cutoff "$OUTDIR/small_vsearch_r1.cutoffs.json"
  log "blast cutoffs:"
  get_cutoff "$OUTDIR/small_blast_r1.cutoffs.json"

  log "vsearch times (s): ${vsearch_times[*]}"
  log "blast   times (s): ${blast_times[*]}"

  local vs_csv1 bl_csv1
  vs_csv1=$(IFS=,; echo "${vsearch_times[*]}")
  bl_csv1=$(IFS=,; echo "${blast_times[*]}")

  python3 -c "
vs = [$vs_csv1]
bl = [$bl_csv1]
vs_mean = sum(vs) / len(vs)
bl_mean = sum(bl) / len(bl)
ratio = bl_mean / vs_mean if vs_mean > 0 else float('inf')
print(f'  vsearch mean: {vs_mean:.2f}s')
print(f'  blast   mean: {bl_mean:.2f}s')
print(f'  ratio (blast/vsearch): {ratio:.2f}x')
"
}

# ── Large test: 620 sequences, local by genus ────────────────────────────────
run_large() {
  bold "\n════ Large dataset: 620 sequences (local by genus) ════"

  local FASTA="tmp/bench/local_many.fasta"
  local CLASS="tmp/bench/local_many.classification"
  local RANK="species"
  local HIGHER_RANK="genus"
  local START_T="0.8"
  local END_T="1.0"

  local vsearch_times=()
  local blast_times=()

  for rep in $(seq 1 "$REPS"); do
    log "Repetition $rep/$REPS"

    rm -f "$OUTDIR"/large_vsearch_r${rep}.*
    rm -f "$OUTDIR"/large_blast_r${rep}.*

    run_timed "large_vsearch_r${rep}" R/predict.R \
      -i "$FASTA" -c "$CLASS" --rank "$RANK" --higher_rank "$HIGHER_RANK" \
      --start_threshold "$START_T" --end_threshold "$END_T" \
      --min_seq_no 2 --min_group_no 2 \
      -o "$OUTDIR" --prefix "large_vsearch_r${rep}" \
      --n_cpus "$N_CPUS" --run_parallel yes --redo yes \
      --tmp_dir "$OUTDIR/tmp_vs"

    run_timed "large_blast_r${rep}" R/predict_blast.R \
      -i "$FASTA" -c "$CLASS" --rank "$RANK" --higher_rank "$HIGHER_RANK" \
      --start_threshold "$START_T" --end_threshold "$END_T" \
      --min_seq_no 2 --min_group_no 2 \
      -o "$OUTDIR" --prefix "large_blast_r${rep}" \
      --n_cpus "$N_CPUS" --run_parallel yes --redo yes \
      --tmp_dir "$OUTDIR/tmp_bl"

    vsearch_times+=("$(cat "$OUTDIR/large_vsearch_r${rep}.time")")
    blast_times+=("$(cat "$OUTDIR/large_blast_r${rep}.time")")
  done

  bold "\n── Large dataset results ──"
  log "vsearch cutoffs (rep 1):"
  get_cutoff "$OUTDIR/large_vsearch_r1.cutoffs.json"
  log "blast cutoffs (rep 1):"
  get_cutoff "$OUTDIR/large_blast_r1.cutoffs.json"

  log "vsearch times (s): ${vsearch_times[*]}"
  log "blast   times (s): ${blast_times[*]}"

  local vs_csv2 bl_csv2
  vs_csv2=$(IFS=,; echo "${vsearch_times[*]}")
  bl_csv2=$(IFS=,; echo "${blast_times[*]}")

  python3 -c "
vs = [$vs_csv2]
bl = [$bl_csv2]
vs_mean = sum(vs) / len(vs)
bl_mean = sum(bl) / len(bl)
ratio = bl_mean / vs_mean if vs_mean > 0 else float('inf')
print(f'  vsearch mean: {vs_mean:.2f}s')
print(f'  blast   mean: {bl_mean:.2f}s')
print(f'  ratio (blast/vsearch): {ratio:.2f}x')
"
}

# ── Main ──────────────────────────────────────────────────────────────────────
bold "╔══════════════════════════════════════════════════════════╗"
bold "║  Speed test: vsearch --allpairs_global vs blastn-short  ║"
bold "╠══════════════════════════════════════════════════════════╣"
log  "CPUs: $N_CPUS   Reps: $REPS   Mode: $MODE"
bold "╚══════════════════════════════════════════════════════════╝"

case "$MODE" in
  small) run_small ;;
  large) run_large ;;
  both)  run_small; run_large ;;
  *)     echo "Usage: $0 [small|large|both]"; exit 1 ;;
esac

bold "\n════ Speed test complete ════"
log "Results in: $OUTDIR/"
