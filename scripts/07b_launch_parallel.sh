#!/bin/bash
# Script name:  07b_launch_parallel.sh
# Description:  Submit 07c_predict_cutoffs_parallel.sh as one SLURM job per ITS
#               region (full_ITS, ITS1, ITS2), so the regions run in parallel.
#
# Usage:        bash scripts/07b_launch_parallel.sh
# Note:         Must be run from the project root directory.

set -euo pipefail

SCRIPT="scripts/07c_predict_cutoffs_parallel.sh"

if [[ ! -f "$SCRIPT" ]]; then
    echo "ERROR: Worker script not found: $SCRIPT" >&2
    exit 1
fi

for region in full_ITS ITS1 ITS2; do
    echo "Submitting ${region}..."
    sbatch --job-name="predict_${region}" "$SCRIPT" "$region"
done

echo ""
echo "All three regions submitted. Monitor with:  squeue -u \$USER"
