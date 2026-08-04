#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Script: prep_featurecounts.annotation.sh
#
# Purpose
# -------
# Prepares a flattened DEXSeq annotation and matching featureCounts GTF
#      using dexseq_prepare_annotation2.py (if outputs do not already exist).
#
# Notes
# -----
# - DEXSeq recommends summarizeOverlaps() in the main vignette, but also points
#   users to the featureCounts-based workflow described by Vivek Bhardwaj.
#
# Inputs
# ------
# resources/genome/
#   Reference GTF annotation used for alignment / counting
#
# Outputs
# -------
# resources/dexseq_annotation/
#   {GTF_BASE}.flat.featurecounts.gff
#   {GTF_BASE}.flat.featurecounts.gtf
#   dexseq_prepare_annotation2.py
#
# Dependencies
# ------------
# External tools:
#   python3
#   curl
#
# Internal utilities:
#   scripts/utils/load_project.sh
#
# ==============================================================================

# -------------------------------
# Helper Functions
# -------------------------------
die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_dir() { [[ -d "$1" ]] || die "Missing directory: $1"; }
require_file() { [[ -s "$1" ]] || die "Missing or empty file: $1"; }

require_tool() {
    local tool="$1"
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "[ERROR] Required tool '$tool' not found in PATH." >&2
        echo "[ERROR] Activate the correct conda environment or install it." >&2
        exit 1
    fi
}

# -------------------------------
# Load Project Constants / Directory Paths
# -------------------------------
# shellcheck disable=SC1090
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../utils/load_project.sh"

: "${GENOME_GTF:?GENOME_GTF not set in environment}"

require_dir "${PROJECT_ROOT}"
require_dir "${RESULTS_DIR}"
require_file "${GENOME_GTF}"

require_tool python3
require_tool curl

# -------------------------------
# Set Paths / Parameters
# -------------------------------
ANNOT_DIR="${PROJECT_ROOT}/resources/dexseq_annotation"
SCRIPTS_DIR="${PROJECT_ROOT}/scripts/utils"

GTF_BASE="$(basename "${GENOME_GTF}")"
GTF_BASE="${GTF_BASE%.gtf}"
GTF_BASE="${GTF_BASE%.gz}"

DEXSEQ_GFF="${ANNOT_DIR}/${GTF_BASE}.flat.featurecounts.gff"
DEXSEQ_FC_GTF="${ANNOT_DIR}/${GTF_BASE}.flat.featurecounts.gtf"

DEXSEQ_PREP_SCRIPT="${SCRIPTS_DIR}/dexseq_prepare_annotation2.py"

mkdir -p "${ANNOT_DIR}" "${SCRIPTS_DIR}"

# -------------------------------
# Prepare DEXSeq / featureCounts annotation
# -------------------------------
if [[ -s "${DEXSEQ_GFF}" && -s "${DEXSEQ_FC_GTF}" ]]; then
    echo "[featureCounts] Found existing flattened annotation."
else
    echo "[featureCounts] Preparing flattened DEXSeq annotation..."

    if [[ ! -s "${DEXSEQ_PREP_SCRIPT}" ]]; then
        echo "[featureCounts] Downloading dexseq_prepare_annotation2.py ..."
        curl -L \
            -o "${DEXSEQ_PREP_SCRIPT}" \
            "https://raw.githubusercontent.com/vivekbhr/Subread_to_DEXSeq/master/dexseq_prepare_annotation2.py"
    fi

    require_file "${DEXSEQ_PREP_SCRIPT}"

    python3 "${DEXSEQ_PREP_SCRIPT}" \
        -f "${DEXSEQ_FC_GTF}" \
        "${GENOME_GTF}" \
        "${DEXSEQ_GFF}"

    require_file "${DEXSEQ_GFF}"
    require_file "${DEXSEQ_FC_GTF}"

    echo "[featureCounts] Flattened annotation complete."
fi

echo "[featureCounts] Prepping annotation complete."
echo "  annotation gff: ${DEXSEQ_GFF}"
echo "  annotation gtf: ${DEXSEQ_FC_GTF}"
