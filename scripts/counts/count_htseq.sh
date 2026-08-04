#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Script: count_htseq.sh
#
# Purpose
# -------
# Generate DEXSeq/HTSeq exon-bin count files from STAR BAMs using the Python
# scripts bundled with the DEXSeq R package:
#   1. dexseq_prepare_annotation.py
#   2. dexseq_count.py
#
# This script:
#   - prepares a flattened DEXSeq GFF (if missing)
#   - iterates over samples.tsv
#   - runs dexseq_count.py on each STAR BAM
#
# Assumptions
# -----------
# - paired-end data
# - BAMs are coordinate-sorted (STAR SortedByCoordinate output)
# - unstranded library
#
# Example Execution
# -----------------
# From notebook or command line:
#
#   "$PROJECT_ROOT/scripts/run" bash \
#       "$PROJECT_ROOT/scripts/counts/count_htseq.sh"

#
# ==============================================================================

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
        exit 1
    fi
}

# shellcheck disable=SC1090
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../utils/load_project.sh"

: "${PROJECT_ROOT:?PROJECT_ROOT not set}"
: "${RESULTS_DIR:?RESULTS_DIR not set}"
: "${SAMPLES_TSV:?SAMPLES_TSV not set}"
: "${GENOME_GTF:?GENOME_GTF not set}"

require_dir "${PROJECT_ROOT}"
require_dir "${RESULTS_DIR}"
require_file "${SAMPLES_TSV}"
require_file "${GENOME_GTF}"

require_tool python3
require_tool Rscript

GTF_BASE="$(basename "${GENOME_GTF}")"
GTF_BASE="${GTF_BASE%.gtf}"
GTF_BASE="${GTF_BASE%.gz}"

DEXSEQ_COUNTS_DIR="${RESULTS_DIR}/counts/htseq"
DEXSEQ_FLAT_GFF="${PROJECT_ROOT}/resources/dexseq_annotation/${GTF_BASE}.flat.htseq.gff"

mkdir -p "${DEXSEQ_COUNTS_DIR}"

# Locate DEXSeq python scripts from the installed R package
DEXSEQ_PYTHON_DIR="$(
    Rscript -e 'cat(system.file("python_scripts", package="DEXSeq", mustWork=TRUE))'
)"

DEXSEQ_PREP="${DEXSEQ_PYTHON_DIR}/dexseq_prepare_annotation.py"
DEXSEQ_COUNT="${DEXSEQ_PYTHON_DIR}/dexseq_count.py"

require_file "${DEXSEQ_PREP}"
require_file "${DEXSEQ_COUNT}"

# Build flattened GFF once
if [[ -s "${DEXSEQ_FLAT_GFF}" ]]; then
    echo "[HTSeq] Found flattened annotation: ${DEXSEQ_FLAT_GFF}"
else
    echo "[HTSeq] Building flattened annotation..."
    python3 "${DEXSEQ_PREP}" "${GENOME_GTF}" "${DEXSEQ_FLAT_GFF}"
    require_file "${DEXSEQ_FLAT_GFF}"
fi

MANIFEST="${DEXSEQ_COUNTS_DIR}/count_manifest.tsv"
echo -e "sample_id\tcondition\tcount_file" >"${MANIFEST}"

tail -n +2 "${SAMPLES_TSV}" |
    while IFS=$'\t' read -r sample cond rep r1_rel r2_rel layout stranded excl; do
        [[ -z "${sample}" ]] && continue
        excl="${excl//$'\r'/}"

        if [[ "${excl}" == "TRUE" ]]; then
            echo "[HTSeq] Skipping excluded sample ${sample}"
            continue
        fi

        bam="${RESULTS_DIR}/counts/star/${cond}/${sample}/${sample}.Aligned.sortedByCoord.out.bam"
        require_file "${bam}"

        out_txt="${DEXSEQ_COUNTS_DIR}/${sample}.dexseq.txt"

        if [[ -s "${out_txt}" ]]; then
            echo "[HTSeq] Existing count file found; skipping ${sample}"
        else
            echo "[HTSeq] Counting ${sample}"
            python3 "${DEXSEQ_COUNT}" \
                -f bam \
                -r pos \
                -p yes \
                -s no \
                "${DEXSEQ_FLAT_GFF}" \
                "${bam}" \
                "${out_txt}"
        fi

        require_file "${out_txt}"
        echo -e "${sample}\t${cond}\t${out_txt}" >>"${MANIFEST}"
    done

echo "[HTSeq] Done."
echo "  flattened gff: ${DEXSEQ_FLAT_GFF}"
echo "  counts dir:    ${DEXSEQ_COUNTS_DIR}"
echo "  manifest:      ${MANIFEST}"
