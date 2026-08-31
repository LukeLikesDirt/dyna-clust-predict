#!/bin/bash
# sync_to_sibling.sh — Materialise this repo's manifest-listed files (as
# committed at HEAD) into a sibling repo's working tree.
#
# Reads tools/sync_manifest.tsv, and for every row where the given target
# column is "yes", writes that file's committed content (git show HEAD:<path>)
# into <target_repo_dir>/<path>, creating parent directories as needed.
#
# Deliberately does NOT git add/commit/push -- it only materialises files.
# The caller decides what to do next: the initial catch-up sync (run by hand,
# reviewed, then committed) and the GitHub Action (which stages, commits,
# and opens a PR) both use this same script as the one place the actual
# file-copying logic lives.
#
# Usage:
#   tools/sync_to_sibling.sh <am|ecm> <target_repo_dir>
#
# Must be run from the dyna-clust-predict project root (reads
# tools/sync_manifest.tsv and git show HEAD:<path> relative to cwd).

set -euo pipefail

target="${1:?Usage: sync_to_sibling.sh <am|ecm> <target_repo_dir>}"
target_dir="${2:?Usage: sync_to_sibling.sh <am|ecm> <target_repo_dir>}"

manifest="tools/sync_manifest.tsv"

if [[ ! -f "$manifest" ]]; then
    echo "ERROR: manifest not found: $manifest (run from the project root)" >&2
    exit 1
fi
if [[ "$target" != "am" && "$target" != "ecm" ]]; then
    echo "ERROR: target must be 'am' or 'ecm', got: $target" >&2
    exit 1
fi
if [[ ! -d "$target_dir/.git" ]]; then
    echo "ERROR: not a git repo: $target_dir" >&2
    exit 1
fi

# Column index: 2 = am, 3 = ecm (column 1 is path)
col=2
[[ "$target" == "ecm" ]] && col=3

n_synced=0
n_skipped=0

while IFS=$'\t' read -r path am_flag ecm_flag; do
    [[ "$path" =~ ^#.*$ || -z "$path" ]] && continue

    flag="$am_flag"
    [[ "$target" == "ecm" ]] && flag="$ecm_flag"

    if [[ "$flag" != "yes" ]]; then
        n_skipped=$((n_skipped + 1))
        continue
    fi

    if ! git cat-file -e "HEAD:$path" 2>/dev/null; then
        echo "WARNING: $path is in the manifest but not tracked at HEAD -- skipping" >&2
        continue
    fi

    dest="$target_dir/$path"
    mkdir -p "$(dirname "$dest")"
    git show "HEAD:$path" > "$dest"
    echo "  synced: $path -> $dest"
    n_synced=$((n_synced + 1))
done < "$manifest"

echo ""
echo "Done. $n_synced file(s) synced to $target ($target_dir), $n_skipped not applicable to $target."
