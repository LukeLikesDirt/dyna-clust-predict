#!/bin/bash
# check_drift.sh — Report how far each sibling repo has drifted from this
# repo's manifest-listed files, without needing to wait on CI. Works offline
# (e.g. on the HPC), as long as the sibling checkouts are available locally.
#
# For every path in tools/sync_manifest.tsv applicable to a given sibling,
# diffs this repo's committed content (git show HEAD:<path>) against that
# path in the sibling's working tree, and reports MISSING / DIFFERS / OK.
#
# Usage:
#   tools/check_drift.sh                 # checks both siblings at their
#                                         # default sibling-directory paths
#   tools/check_drift.sh am              # checks only dyna-clust-predict-am
#   tools/check_drift.sh ecm             # checks only dyna-clust-predict-ecm
#
# Override the sibling checkout locations (default: ../dyna-clust-predict-am
# and ../dyna-clust-predict-ecm, i.e. siblings of this repo's own directory)
# via environment variables:
#   DYNA_CLUST_AM_DIR=/path/to/am tools/check_drift.sh
#   DYNA_CLUST_ECM_DIR=/path/to/ecm tools/check_drift.sh
#
# Must be run from the dyna-clust-predict project root.

set -uo pipefail

manifest="tools/sync_manifest.tsv"
if [[ ! -f "$manifest" ]]; then
    echo "ERROR: manifest not found: $manifest (run from the project root)" >&2
    exit 1
fi

am_dir="${DYNA_CLUST_AM_DIR:-../dyna-clust-predict-am}"
ecm_dir="${DYNA_CLUST_ECM_DIR:-../dyna-clust-predict-ecm}"

targets=("am" "ecm")
if [[ $# -ge 1 ]]; then
    case "$1" in
        am|ecm) targets=("$1") ;;
        *) echo "ERROR: argument must be 'am' or 'ecm'" >&2; exit 1 ;;
    esac
fi

overall_exit=0

check_target() {
    local target="$1" target_dir="$2"

    echo "=== $target  ($target_dir) ==="
    if [[ ! -d "$target_dir" ]]; then
        local target_upper
        target_upper=$(echo "$target" | tr '[:lower:]' '[:upper:]')
        echo "  SKIPPED: directory not found (set DYNA_CLUST_${target_upper}_DIR to override)"
        echo ""
        return
    fi

    local n_ok=0 n_diff=0 n_missing=0

    while IFS=$'\t' read -r path am_flag ecm_flag; do
        [[ "$path" =~ ^#.*$ || -z "$path" ]] && continue

        local flag="$am_flag"
        [[ "$target" == "ecm" ]] && flag="$ecm_flag"
        [[ "$flag" != "yes" ]] && continue

        local sibling_file="$target_dir/$path"
        if [[ ! -f "$sibling_file" ]]; then
            echo "  MISSING   $path"
            n_missing=$((n_missing + 1))
            overall_exit=1
            continue
        fi

        local n_changed
        n_changed=$(git diff --no-index --numstat <(git show "HEAD:$path") "$sibling_file" 2>/dev/null \
                    | awk '{print $1+$2}')

        if [[ -z "$n_changed" || "$n_changed" == "0" ]]; then
            n_ok=$((n_ok + 1))
        else
            echo "  DIFFERS   $path  ($n_changed line(s) different)"
            n_diff=$((n_diff + 1))
            overall_exit=1
        fi
    done < "$manifest"

    echo "  -> $n_ok in sync, $n_diff differ, $n_missing missing"
    echo ""
}

for t in "${targets[@]}"; do
    if [[ "$t" == "am" ]]; then check_target "am" "$am_dir"; fi
    if [[ "$t" == "ecm" ]]; then check_target "ecm" "$ecm_dir"; fi
done

exit "$overall_exit"
