#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=0-01:00:00
#SBATCH --partition=short
#SBATCH --output=logs/%x.%j.out

# Script name:  01_reformat_ITS.sh
# Description:  Download the EUKARYOME ITS v2.0 database, reformat sequence
#               headers with reformat.R, extract the taxonomy classification
#               table, and canonicalise it so every named taxon has exactly
#               one parent (canonicalise_lineage.R).
# Note:         This script must be run from the project root directory.

# =============================================================================
# PARAMETER SETUP
# =============================================================================

# EUKARYOME download URL (EUKARYOME General EUK ITS v2.0)
readonly EUKARYOME_URL="https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_ITS_v2.0.zip"
readonly DOWNLOAD_FILE="./tmp/General_EUK_ITS_v2.0.zip"
readonly EXTRACTED_FASTA="./tmp/General_EUK_ITS_v2.0.fasta"

# Output files
readonly OUT_FASTA="./tmp/eukaryome_ITS.fasta"
readonly OUT_CLASS="./tmp/eukaryome_ITS.classification"

# Helper R scripts
readonly REFORMAT="./R/reformat.R"
readonly CANONICALISE="./R/canonicalise_lineage.R"

# Labels that had more than one parent in EUKARYOME, and the parent chosen
readonly LINEAGE_MANIFEST="./data/lineage_canonicalisation.txt"

# =============================================================================
# DIRECTORY SETUP
# =============================================================================

mkdir -p ./data/full_ITS tmp

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

echo "Activating conda environment..."
source ~/.bashrc
conda activate dyna_clust_predict

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if [[ ! -f "$REFORMAT" ]]; then
    echo "ERROR: R script not found: $REFORMAT" >&2
    exit 1
fi

if [[ ! -f "$CANONICALISE" ]]; then
    echo "ERROR: R script not found: $CANONICALISE" >&2
    exit 1
fi

# =============================================================================
# FILE DOWNLOAD: EUKARYOME
# =============================================================================

echo ""
echo "=== DOWNLOADING EUKARYOME DATABASE ==="
echo "$(date)"
echo "URL: $EUKARYOME_URL"

if ! curl -L -o "$DOWNLOAD_FILE" "$EUKARYOME_URL"; then
    echo "ERROR: Failed to download file from $EUKARYOME_URL" >&2
    exit 1
fi

echo "Unzipping downloaded file..."
if ! 7z x "$DOWNLOAD_FILE" -o"./tmp/" -y; then
    echo "ERROR: Failed to unzip $DOWNLOAD_FILE" >&2
    exit 1
fi

echo "Download and extraction completed successfully."
echo ""

# =============================================================================
# REFORMAT HEADERS
# =============================================================================

echo "=== REFORMATTING HEADERS ==="
echo "$(date)"

if [[ ! -f "$EXTRACTED_FASTA" ]]; then
    echo "ERROR: Extracted FASTA not found: $EXTRACTED_FASTA" >&2
    exit 1
fi

echo "Running reformat.R..."
Rscript "$REFORMAT" \
    --fasta_in        "$EXTRACTED_FASTA" \
    --fasta_out       "$OUT_FASTA" \
    --classification_out "$OUT_CLASS"

if [[ $? -ne 0 ]]; then
    echo "ERROR: reformat.R failed." >&2
    exit 1
fi

echo ""
echo "Reformatted FASTA written to      : $OUT_FASTA"
echo "Classification table written to   : $OUT_CLASS"
echo ""

# =============================================================================
# CANONICALISE THE LINEAGE
# =============================================================================
# EUKARYOME places some names under more than one parent (e.g. Cladosporiaceae
# under Cladosporiales and Dothideales). Every named label is given its modal
# identified parent, so the classification is a tree before dereplication, and
# 08_consolidate_cutoffs.sh no longer has to guess a lineage ("inconsistent
# parent ... modal parent used for lineage lookup").

echo "=== CANONICALISING THE LINEAGE ==="
echo "$(date)"

Rscript "$CANONICALISE" \
    --classification_in  "$OUT_CLASS" \
    --classification_out "$OUT_CLASS" \
    --manifest_out       "$LINEAGE_MANIFEST"

if [[ $? -ne 0 ]]; then
    echo "ERROR: canonicalise_lineage.R failed." >&2
    exit 1
fi

echo ""
echo "Canonical classification written to: $OUT_CLASS"
echo "Re-parented labels listed in       : $LINEAGE_MANIFEST"
echo ""

# =============================================================================
# CLEANUP
# =============================================================================

echo "=== CLEANUP ==="
echo "Removing downloaded and extracted files from tmp/..."
rm -f "$DOWNLOAD_FILE" "$EXTRACTED_FASTA"

echo ""
echo "=== PIPELINE COMPLETED SUCCESSFULLY ==="
echo "$(date)"

conda deactivate
