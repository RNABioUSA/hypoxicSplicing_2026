#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Script: count_salmon.sh
#
# Purpose
# -------
# Perform transcript quantification using Salmon for all samples listed
# in the project sample sheet. The script:
#
#   1. Builds a decoy-aware Salmon index from the reference genome and GTF
#      (if the index does not already exist).
#   2. Iterates through the sample metadata table.
#   3. Quantifies transcript abundance for each sample.
#   4. Generates a coldata file mapping sample names to quant.sf files
#      for downstream analysis.
#
# The script avoids recomputing results when outputs already exist.
#
# Workflow Overview
# -----------------
# FASTQ reads
#      │
#      ▼
# salmon index
#      │
#      ▼
# salmon quant
#      │
#      ▼
# quant.sf
#
# Inputs
# ------
# reads/
#   Raw FASTQ files referenced in samples.tsv
#
# config/metadata/samples.tsv
#   Sample metadata table containing:
#
#       sample_id
#       condition
#       replicate
#       read1_rel
#       read2_rel
#       layout
#       stranded
#       exclude
#
# resources/genome/
#   Reference genome FASTA and GTF annotation
#
# Outputs
# -------
# resources/salmon_index/
#   Decoy-aware Salmon index
#
# results/counts/salmon/
#   <condition>/<sample>/
#       quant.sf
#       cmd_info.json
#
#   coldata.tsv
#
# Dependencies
# ------------
# External tools (expected in PATH via conda environment):
#
#   salmon
#   gffread
#
# Internal utilities:
#
#   scripts/utils/load_project.sh
#
# Environment Variables
# ---------------------
# Set via config/environment/analysis_environment.env
#
#   PROJECT_ROOT
#   THREADS
#   GENOME_BUILD
#   GENOME_FA
#   GENOME_GTF
#
# Example Execution
# -----------------
# From notebook or command line:
#
#   "$PROJECT_ROOT/scripts/run" bash \
#       "$PROJECT_ROOT/scripts/counts/count_salmon.sh"
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

: "${THREADS:?THREADS not set in environment}"
: "${GENOME_BUILD:?GENOME_BUILD not set in environment}"
: "${GENOME_FA:?GENOME_FA not set in environment}"
: "${GENOME_GTF:?GENOME_GTF not set in environment}"

require_dir "${PROJECT_ROOT}"
require_dir "${READS_ROOT}"
require_dir "${RESULTS_DIR}"
require_file "${SAMPLES_TSV}"

require_tool salmon
require_tool gffread

SALMON_INDEX="${PROJECT_ROOT}/resources/salmon_index"
INDEX_SENTINEL="${SALMON_INDEX}/versionInfo.json"

OUTDIR="${RESULTS_DIR}/counts/salmon"
COLDATA_OUT="${OUTDIR}/coldata.tsv"

mkdir -p "${OUTDIR}"

DONE_FILE="${OUTDIR}/.salmon_quant.done"

# -------------------------------
# Global Skip Check
# -------------------------------
if [[ -s "${DONE_FILE}" ]]; then
    echo "[salmon] Quantification already complete (${DONE_FILE}); skipping."
    exit 0
fi

# -------------------------------
# Set Parameters
# -------------------------------
SALMON_KMER=31
SALMON_OPTS="--libType A --rangeFactorizationBins 4 --gcBias --seqBias --recoverOrphans --validateMappings"
read -r -a SALMON_OPTS_ARR <<<"${SALMON_OPTS}"

# -------------------------------
# Build Index (skip if exists)
# -------------------------------
if [[ ! -s "${INDEX_SENTINEL}" ]]; then
    echo "[salmon] Index not found; building decoy-aware gentrome + index..."

    require_dir "${GENOME_DIR}"
    require_file "${GENOME_FA}"
    require_file "${GENOME_GTF}"

    mkdir -p "${SALMON_INDEX}"

    # Derived names based on GTF (captures version info)
    # Example: gencode.v45.primary_assembly.basic.annotation.gtf -> gencode.v45.primary_assembly.basic.annotation
    gtf_base="$(basename "${GENOME_GTF}" .gtf)"
    gtf_base="${gtf_base%.gz}"

    TRANSCRIPTS_FA="${GENOME_DIR}/${gtf_base}.transcripts.fa"
    TRANSCRIPTS_FA_GZ="${TRANSCRIPTS_FA}.gz"
    DECOYS="${GENOME_DIR}/${GENOME_BUILD}.salmon_decoys.txt"
    GENTROME="${GENOME_DIR}/${gtf_base}.${GENOME_BUILD}.gentrome.fa.gz"

    GENOME_FA_FOR_GFFREAD="${GENOME_FA}"
    TMP_GENOME=""

    if [[ "${GENOME_FA}" == *.gz ]]; then
        TMP_GENOME="$(mktemp --suffix=.fa)"
        echo "[salmon] Decompressing genome FASTA for gffread -> ${TMP_GENOME}"
        gunzip -c "${GENOME_FA}" >"${TMP_GENOME}"
        GENOME_FA_FOR_GFFREAD="${TMP_GENOME}"
    fi

    echo "[salmon] Generating transcripts FASTA with gffread..."
    gffread -w "${TRANSCRIPTS_FA}" -g "${GENOME_FA_FOR_GFFREAD}" "${GENOME_GTF}"
    gzip -c "${TRANSCRIPTS_FA}" >"${TRANSCRIPTS_FA_GZ}"

    echo "[salmon] Generating decoys list..."
    if [[ "${GENOME_FA}" == *.gz ]]; then
        grep "^>" <(gunzip -c "${GENOME_FA}") | cut -d " " -f 1 | sed 's/>//g' >"${DECOYS}"
    else
        grep "^>" "${GENOME_FA}" | cut -d " " -f 1 | sed 's/>//g' >"${DECOYS}"
    fi

    echo "[salmon] Building gentrome..."
    if [[ "${GENOME_FA}" == *.gz ]]; then
        cat "${TRANSCRIPTS_FA_GZ}" "${GENOME_FA}" >"${GENTROME}"
    else
        TMP_GENOME_GZ="$(mktemp --suffix=.fa.gz)"
        gzip -c "${GENOME_FA}" >"${TMP_GENOME_GZ}"
        cat "${TRANSCRIPTS_FA_GZ}" "${TMP_GENOME_GZ}" >"${GENTROME}"
        rm -f "${TMP_GENOME_GZ}"
    fi

    [[ -n "${TMP_GENOME}" ]] && rm -f "${TMP_GENOME}"

    echo "[salmon] Building index at ${SALMON_INDEX} ..."
    salmon index \
        -p "${THREADS}" \
        --gencode \
        --keepDuplicates \
        -k "${SALMON_KMER}" \
        -t "${GENTROME}" \
        -d "${DECOYS}" \
        -i "${SALMON_INDEX}"

    echo "[salmon] Index build complete."
else
    echo "[salmon] Found existing index: ${SALMON_INDEX}"
fi

# -------------------------------
# Quantify per sample from samples.tsv
# Required columns: sample_id, condition, read1_rel, read2_rel
# Optional: exclude (TRUE/FALSE)
# -------------------------------
echo -e "names\tfiles\tcondition" >"${COLDATA_OUT}"

tail -n +2 "${SAMPLES_TSV}" |
    while IFS=$'\t' read -r sample cond rep r1_rel r2_rel layout stranded excl; do

        [[ -z "${sample}" ]] && continue
        excl="${excl//$'\r'/}"

        if [[ "${excl}" == "TRUE" ]]; then
            echo "[salmon] Skipping excluded sample ${sample}"
            continue
        fi

        R1="${READS_ROOT}/${r1_rel}"
        R2="${READS_ROOT}/${r2_rel}"
        require_file "${R1}"
        require_file "${R2}"

        sample_out="${OUTDIR}/${cond}/${sample}"
        sample_out_rel="results/counts/salmon/${cond}/${sample}"
        mkdir -p "${sample_out}"

        if [[ -s "${sample_out}/quant.sf" && -s "${sample_out}/cmd_info.json" ]]; then
            echo "[salmon] Existing quant found; skipping ${sample}"
        else
            echo "[salmon] Quantifying ${sample} (${cond})"
            salmon quant \
                "${SALMON_OPTS_ARR[@]}" \
                -p "${THREADS}" \
                -i "${SALMON_INDEX}" \
                -1 "${R1}" -2 "${R2}" \
                -o "${sample_out}"
        fi

        echo -e "${sample}\t${sample_out_rel}/quant.sf\t${cond}" >>"${COLDATA_OUT}"
    done

date -Is >"${DONE_FILE}"
echo "[salmon] Done."
echo "  quants:     ${OUTDIR}/<sample_id>/quant.sf"
echo "  coldata:    ${COLDATA_OUT}"
