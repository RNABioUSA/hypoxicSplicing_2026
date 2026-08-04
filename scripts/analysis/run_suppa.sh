#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Script: run_suppa.sh
#
# Purpose
# -------
# Run SUPPA2 alternative splicing analysis using transcript TPM estimates.
#
# This script supports two execution modes:
#
#   1. Standard mode
#      -------------
#      Builds isoform TPM matrices directly from Salmon quantification outputs,
#      splits them by condition, calculates PSI values, and runs SUPPA2
#      diffSplice for each comparison against the reference condition.
#
#      Output directory:
#
#          results/analysis/suppa/standard/<reference_condition>/
#
#   2. Harmonized-input mode
#      ---------------------
#      Uses pre-filtered isoform TPM files supplied by the caller, for example
#      TPM files filtered to a DRIMSeq-derived harmonized transcript universe.
#
#      In this mode, the script skips:
#
#          - Salmon TPM matrix construction
#          - TPM splitting
#          - generic per-condition PSI generation
#
#      Instead, contrast-specific filtered TPM files are copied into the SUPPA
#      output directory and PSI values are computed per contrast immediately
#      before diffSplice.
#
#      This is important when the transcript/gene universe differs by contrast.
#
#      Required environment variables:
#
#          SUPPA_RUN_TAG=harmonized
#          SUPPA_INPUT_ISO_TPM_DIR=/path/to/filtered/iso_tpm/files
#
#      Expected harmonized input filenames:
#
#          iso_tpm_<REF>v<COND>_filtered.txt
#          iso_tpm_<COND>_filtered.txt
#
#      Example:
#
#          iso_tpm_C1vH1_filtered.txt
#          iso_tpm_H1_filtered.txt
#
#      Output directory:
#
#          results/analysis/suppa/harmonized/<reference_condition>/
#
# Workflow Overview: Standard Mode
# --------------------------------
# Salmon quant.sf outputs
#        │
#        ▼
# multipleFieldSelection.py
#        │
#        ▼
# iso_tpm.txt
#        │
#        ▼
# split_iso_tpm
#        │
#        ▼
# condition TPM matrices
#        │
#        ▼
# suppa.py psiPerEvent
#        │
#        ▼
# PSI matrices
#        │
#        ▼
# suppa.py diffSplice
#
# Workflow Overview: Harmonized Mode
# ----------------------------------
# filtered contrast-specific isoTPM files
#        │
#        ▼
# copy into harmonized SUPPA output directory
#        │
#        ▼
# suppa.py psiPerEvent per contrast
#        │
#        ▼
# suppa.py diffSplice
#
# Inputs
# ------
# resources/genome/
#   Reference genome annotation, defined by GENOME_GTF.
#
# results/counts/salmon/
#   Salmon quantification results, used only in standard mode.
#
# config/metadata/samples.tsv
#   Sample metadata table with columns:
#       sample_id, condition, ..., exclude
#
# Optional harmonized input directory:
#   SUPPA_INPUT_ISO_TPM_DIR
#
# Outputs
# -------
# resources/suppa_events/
#   <GENOME_BUILD>.events.ioe
#   <GENOME_BUILD>.events.ioi
#
# results/analysis/suppa/<run_tag>/<reference_condition>/
#   iso_tpm/
#   events/
#   diffSplice/
#
# Dependencies
# ------------
# External tools expected in PATH:
#
#   suppa.py
#   multipleFieldSelection.py
#   awk
#   python3
#
# Internal utilities:
#
#   scripts/utils/load_project.sh
#
# Environment Variables
# ---------------------
# Required via project environment:
#
#   PROJECT_ROOT
#   RESULTS_DIR
#   GENOME_BUILD
#   GENOME_GTF
#   SAMPLES_TSV
#
# Optional:
#
#   SUPPA_REF_CONDITION      Default: C1
#   SUPPA_RUN_TAG            Default: standard
#   SUPPA_INPUT_ISO_TPM_DIR  Default: empty
#
# Example Standard Execution
# --------------------------
#   "$PROJECT_ROOT/scripts/run" bash \
#       "$PROJECT_ROOT/scripts/analysis/run_suppa.sh" C1
#
# Example Harmonized Execution
# ----------------------------
#   SUPPA_RUN_TAG=harmonized \
#   SUPPA_INPUT_ISO_TPM_DIR="$RESULTS_DIR/analysis/suppa_harmonized_inputs/iso_tpm" \
#   "$PROJECT_ROOT/scripts/run" bash \
#       "$PROJECT_ROOT/scripts/analysis/run_suppa.sh" C1
#
# NOTE:
# SUPPA multipleFieldSelection.py writes a non-standard matrix where the header
# contains sample names only and each data row begins with a transcript ID. The
# custom splitter below preserves this SUPPA-compatible layout.
#
# ==============================================================================


# -------------------------------
# Helper Functions
# -------------------------------
die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_dir() {
    [[ -d "$1" ]] || die "Missing directory: $1"
}

require_file() {
    [[ -s "$1" ]] || die "Missing or empty file: $1"
}

require_tool() {
    local tool="$1"

    if ! command -v "${tool}" >/dev/null 2>&1; then
        echo "[ERROR] Required tool '${tool}' not found in PATH." >&2
        echo "[ERROR] Activate the correct conda environment or install it." >&2
        exit 1
    fi
}

get_samples_csv() {
    local cond="$1"

    awk -F '\t' -v target="${cond}" '
        NR == 1 { next }
        {
            excl = $8
            gsub(/\r/, "", excl)
            if ($2 == target && excl != "TRUE") print $1
        }
    ' "${SAMPLES_TSV}" | paste -sd,
}

split_iso_tpm() {
    local infile="$1"
    local samples1_csv="$2"
    local samples2_csv="$3"
    local out1="$4"
    local out2="$5"

    python3 - "$infile" "$samples1_csv" "$samples2_csv" "$out1" "$out2" <<'PY'
import sys

infile, samples1_csv, samples2_csv, out1, out2 = sys.argv[1:]

samples1 = [x for x in samples1_csv.split(",") if x]
samples2 = [x for x in samples2_csv.split(",") if x]

with open(infile, "r", newline="") as fh:
    # Master SUPPA iso_tpm.txt format:
    # header = sample names only
    # data   = transcript_id + values
    header = fh.readline().rstrip("\r\n").split("\t")

    sample_to_data_idx = {name: i + 1 for i, name in enumerate(header)}

    missing1 = [s for s in samples1 if s not in sample_to_data_idx]
    missing2 = [s for s in samples2 if s not in sample_to_data_idx]

    if missing1:
        raise SystemExit(
            f"[split_iso_tpm] Missing samples for condition 1: {', '.join(missing1)}"
        )

    if missing2:
        raise SystemExit(
            f"[split_iso_tpm] Missing samples for condition 2: {', '.join(missing2)}"
        )

    idx1 = [sample_to_data_idx[s] for s in samples1]
    idx2 = [sample_to_data_idx[s] for s in samples2]

    expected_n = len(header) + 1  # transcript_id + all sample columns

    with open(out1, "w", newline="") as f1, open(out2, "w", newline="") as f2:
        # Write header in the same odd SUPPA format: sample names only
        f1.write("\t".join(samples1) + "\n")
        f2.write("\t".join(samples2) + "\n")

        for lineno, line in enumerate(fh, start=2):
            row = line.rstrip("\r\n").split("\t")

            if len(row) != expected_n:
                raise SystemExit(
                    f"[split_iso_tpm] Line {lineno} in {infile} has "
                    f"{len(row)} fields; expected {expected_n}"
                )

            f1.write("\t".join([row[0]] + [row[i] for i in idx1]) + "\n")
            f2.write("\t".join([row[0]] + [row[i] for i in idx2]) + "\n")
PY
}


# -------------------------------
# Load Project Constants / Directory Paths
# -------------------------------
# shellcheck disable=SC1090
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../utils/load_project.sh"

: "${GENOME_BUILD:?GENOME_BUILD not set in environment}"
: "${GENOME_GTF:?GENOME_GTF not set in environment}"

require_dir "${PROJECT_ROOT}"
require_dir "${RESULTS_DIR}"
require_dir "${GENOME_DIR}"
require_file "${GENOME_GTF}"
require_file "${SAMPLES_TSV}"

require_tool suppa.py
require_tool multipleFieldSelection.py
require_tool awk
require_tool python3


# -------------------------------
# Paths / Run Mode
# -------------------------------
SALMON_DIR="${RESULTS_DIR}/counts/salmon"

SUPPA_RESOURCE_DIR="${PROJECT_ROOT}/resources/suppa_events"
EVENTS_IOE_DIR="${SUPPA_RESOURCE_DIR}/events_ioe"
IOE_MASTER="${SUPPA_RESOURCE_DIR}/${GENOME_BUILD}.events.ioe"
IOI_PREFIX="${SUPPA_RESOURCE_DIR}/${GENOME_BUILD}.events"

# -------------------------------
# Parse CLI Arguments
# -------------------------------
SUPPA_REF_CONDITION="${1:-C1}"
SUPPA_RUN_TAG="${2:-standard}"
SUPPA_INPUT_ISO_TPM_DIR="${3:-}"

# Validate arguments
if [[ "${SUPPA_RUN_TAG}" != "standard" && "${SUPPA_RUN_TAG}" != "harmonized" ]]; then
    die "SUPPA_RUN_TAG must be 'standard' or 'harmonized'"
fi

if [[ "${SUPPA_RUN_TAG}" == "harmonized" && -z "${SUPPA_INPUT_ISO_TPM_DIR}" ]]; then
    die "Harmonized mode requires isoTPM input directory as third argument"
fi

if [[ "${SUPPA_RUN_TAG}" == "standard" ]]; then
    OUTDIR="${RESULTS_DIR}/analysis/suppa/standard/${SUPPA_REF_CONDITION}"
else
    OUTDIR="${RESULTS_DIR}/analysis/suppa/${SUPPA_RUN_TAG}/${SUPPA_REF_CONDITION}"
fi

echo "[SUPPA2] Mode: ${SUPPA_RUN_TAG}"

if [[ "${SUPPA_RUN_TAG}" == "harmonized" ]]; then
    echo "[SUPPA2] Using isoTPM input: ${SUPPA_INPUT_ISO_TPM_DIR}"
fi

ISO_TPM_DIR="${OUTDIR}/iso_tpm"
EVENTS_DIR="${OUTDIR}/events"
DIFF_DIR="${OUTDIR}/diffSplice"

DONE_FILE="${OUTDIR}/.suppa_diffsplice.done"

mkdir -p "${SUPPA_RESOURCE_DIR}" "${EVENTS_IOE_DIR}"
mkdir -p "${OUTDIR}" "${ISO_TPM_DIR}" "${EVENTS_DIR}" "${DIFF_DIR}"


# -------------------------------
# Global Skip Check
# -------------------------------
if [[ -s "${DONE_FILE}" ]]; then
    echo "[SUPPA2] Differential splicing already complete (${DONE_FILE}); skipping."
    exit 0
fi


# -------------------------------
# Validate Core Inputs
# -------------------------------
require_dir "${SALMON_DIR}"

if [[ -n "${SUPPA_INPUT_ISO_TPM_DIR}" ]]; then
    require_dir "${SUPPA_INPUT_ISO_TPM_DIR}"
fi

echo "[SUPPA2] PROJECT_ROOT=${PROJECT_ROOT}"
echo "[SUPPA2] RESULTS_DIR=${RESULTS_DIR}"
echo "[SUPPA2] SALMON_DIR=${SALMON_DIR}"
echo "[SUPPA2] OUTDIR=${OUTDIR}"
echo "[SUPPA2] Reference condition=${SUPPA_REF_CONDITION}"
echo "[SUPPA2] Run tag=${SUPPA_RUN_TAG}"

if [[ -n "${SUPPA_INPUT_ISO_TPM_DIR}" ]]; then
    echo "[SUPPA2] Harmonized isoTPM input dir=${SUPPA_INPUT_ISO_TPM_DIR}"
fi


# -------------------------------
# Generate SUPPA Event Annotations
# -------------------------------
if [[ -s "${IOE_MASTER}" ]]; then
    echo "[SUPPA2] Found existing event annotation: ${IOE_MASTER}"
else
    echo "[SUPPA2] Event annotation not found; generating from ${GENOME_GTF}"

    suppa.py generateEvents \
        -i "${GENOME_GTF}" \
        -o "${EVENTS_IOE_DIR}/events" \
        -f ioe \
        -e SE SS MX RI FL

    suppa.py generateEvents \
        -i "${GENOME_GTF}" \
        -o "${IOI_PREFIX}" \
        -f ioi

    ioe_files=("${EVENTS_IOE_DIR}"/*.ioe)
    [[ -e "${ioe_files[0]}" ]] || die "No .ioe files generated in ${EVENTS_IOE_DIR}"

    echo "[SUPPA2] Merging per-event IOE files -> ${IOE_MASTER}"

    awk '
        FNR == 1 && NR != 1 { next }
        { print }
    ' "${ioe_files[@]}" >"${IOE_MASTER}"

    require_file "${IOE_MASTER}"
fi


# -------------------------------
# Parse Samples
# Assumes:
#   column 1 = sample_id
#   column 2 = condition
#   column 8 = exclude TRUE/FALSE
# -------------------------------
mapfile -t CONDITIONS < <(
    awk -F '\t' '
        NR == 1 { next }
        {
            excl = $8
            gsub(/\r/, "", excl)
            if (excl != "TRUE" && $2 != "") print $2
        }
    ' "${SAMPLES_TSV}" | sort -u
)

((${#CONDITIONS[@]} > 0)) || die "No non-excluded conditions found in ${SAMPLES_TSV}"

printf '%s\n' "${CONDITIONS[@]}" | grep -qx "${SUPPA_REF_CONDITION}" ||
    die "Reference condition '${SUPPA_REF_CONDITION}' not found in ${SAMPLES_TSV}"

REF_SAMPLE_CSV="$(get_samples_csv "${SUPPA_REF_CONDITION}")"
[[ -n "${REF_SAMPLE_CSV}" ]] || die "No non-excluded samples found for ${SUPPA_REF_CONDITION}"


# -------------------------------
# Standard Mode:
# Build Master TPM Matrix and Split by Condition
# -------------------------------
if [[ -z "${SUPPA_INPUT_ISO_TPM_DIR}" ]]; then
    ISO_TPM_ALL="${ISO_TPM_DIR}/iso_tpm.txt"
    REF_ISO="${ISO_TPM_DIR}/iso_tpm_${SUPPA_REF_CONDITION}.txt"
    REF_PSI="${EVENTS_DIR}/events_${SUPPA_REF_CONDITION}.psi"

    if [[ -s "${ISO_TPM_ALL}" ]]; then
        echo "[SUPPA2] Found existing master isoform TPM matrix: ${ISO_TPM_ALL}"
    else
        echo "[SUPPA2] Building master isoform TPM matrix from Salmon outputs..."

        quant_files=("${SALMON_DIR}"/*/*/quant.sf)
        [[ -e "${quant_files[0]}" ]] || die "No Salmon quant.sf files found under ${SALMON_DIR}"

        multipleFieldSelection.py \
            -i "${quant_files[@]}" \
            -k 1 \
            -f 4 \
            -o "${ISO_TPM_ALL}"

        require_file "${ISO_TPM_ALL}"
    fi

    for cond in "${CONDITIONS[@]}"; do
        [[ "${cond}" == "${SUPPA_REF_CONDITION}" ]] && continue

        COND_SAMPLE_CSV="$(get_samples_csv "${cond}")"
        [[ -n "${COND_SAMPLE_CSV}" ]] || continue

        COND_ISO="${ISO_TPM_DIR}/iso_tpm_${cond}.txt"

        if [[ -s "${REF_ISO}" && -s "${COND_ISO}" ]]; then
            echo "[SUPPA2] Existing condition TPM files found; skipping split for ${SUPPA_REF_CONDITION} vs ${cond}"
        else
            echo "[SUPPA2] Splitting TPM matrix for ${SUPPA_REF_CONDITION} vs ${cond}"

            split_iso_tpm \
                "${ISO_TPM_ALL}" \
                "${REF_SAMPLE_CSV}" \
                "${COND_SAMPLE_CSV}" \
                "${REF_ISO}" \
                "${COND_ISO}"
        fi

        require_file "${REF_ISO}"
        require_file "${COND_ISO}"
    done
else
    echo "[SUPPA2] Harmonized isoTPM input provided; skipping Salmon TPM build and split."
fi


# -------------------------------
# Standard Mode:
# Compute PSI Per Event for Each Condition
# -------------------------------
if [[ -z "${SUPPA_INPUT_ISO_TPM_DIR}" ]]; then
    for cond in "${CONDITIONS[@]}"; do
        ISO_FILE="${ISO_TPM_DIR}/iso_tpm_${cond}.txt"
        PSI_PREFIX="${EVENTS_DIR}/events_${cond}"
        PSI_FILE="${PSI_PREFIX}.psi"

        [[ -s "${ISO_FILE}" ]] || continue

        if [[ -s "${PSI_FILE}" ]]; then
            echo "[SUPPA2] Existing PSI file found; skipping ${cond}"
            continue
        fi

        echo "[SUPPA2] Inspecting ${ISO_FILE}"
        awk -F '\t' '
            NR == 1 {
                print "[suppa] header NF=" NF, "|" $0 "|"
            }
            NR == 118012 {
                print "[suppa] target NF=" NF, "|" $0 "|"
            }
        ' "${ISO_FILE}"

        echo "[SUPPA2] Calculating PSI for ${cond}"

        suppa.py psiPerEvent \
            -i "${IOE_MASTER}" \
            -e "${ISO_FILE}" \
            -o "${PSI_PREFIX}"

        require_file "${PSI_FILE}"
    done
else
    echo "[SUPPA2] Harmonized isoTPM input provided; PSI will be computed per contrast."
fi


# -------------------------------
# Run diffSplice
#
# Standard mode:
#   Uses one reference isoTPM and one condition isoTPM per condition.
#
# Harmonized mode:
#   Uses contrast-specific filenames to avoid accidental reuse across contrasts:
#
#       iso_tpm_C1vH1.txt
#       iso_tpm_H1vC1.txt
#
#   PSI files are also contrast-specific:
#
#       events_C1vH1.psi
#       events_H1vC1.psi
# -------------------------------
for cond in "${CONDITIONS[@]}"; do
    [[ "${cond}" == "${SUPPA_REF_CONDITION}" ]] && continue

    OUT_PREFIX="${DIFF_DIR}/diffSplice_${cond}"

    if ls "${OUT_PREFIX}"* >/dev/null 2>&1; then
        echo "[SUPPA2] Existing diffSplice outputs found; skipping ${SUPPA_REF_CONDITION} vs ${cond}"
        continue
    fi

    if [[ -n "${SUPPA_INPUT_ISO_TPM_DIR}" ]]; then
        SRC_REF_ISO="${SUPPA_INPUT_ISO_TPM_DIR}/iso_tpm_${SUPPA_REF_CONDITION}v${cond}_filtered.txt"
        SRC_COND_ISO="${SUPPA_INPUT_ISO_TPM_DIR}/iso_tpm_${cond}_filtered.txt"

        require_file "${SRC_REF_ISO}"
        require_file "${SRC_COND_ISO}"

        REF_ISO="${ISO_TPM_DIR}/iso_tpm_${SUPPA_REF_CONDITION}v${cond}.txt"
        COND_ISO="${ISO_TPM_DIR}/iso_tpm_${cond}v${SUPPA_REF_CONDITION}.txt"

        if [[ ! -s "${REF_ISO}" ]]; then
            echo "[SUPPA2] Copying harmonized reference isoTPM: ${SRC_REF_ISO}"
            cp "${SRC_REF_ISO}" "${REF_ISO}"
        else
            echo "[SUPPA2] Existing harmonized reference isoTPM found; keeping ${REF_ISO}"
        fi

        if [[ ! -s "${COND_ISO}" ]]; then
            echo "[SUPPA2] Copying harmonized condition isoTPM: ${SRC_COND_ISO}"
            cp "${SRC_COND_ISO}" "${COND_ISO}"
        else
            echo "[SUPPA2] Existing harmonized condition isoTPM found; keeping ${COND_ISO}"
        fi

        REF_PSI_PREFIX="${EVENTS_DIR}/events_${SUPPA_REF_CONDITION}v${cond}"
        REF_PSI="${REF_PSI_PREFIX}.psi"

        COND_PSI_PREFIX="${EVENTS_DIR}/events_${cond}v${SUPPA_REF_CONDITION}"
        COND_PSI="${COND_PSI_PREFIX}.psi"
    else
        REF_ISO="${ISO_TPM_DIR}/iso_tpm_${SUPPA_REF_CONDITION}.txt"
        COND_ISO="${ISO_TPM_DIR}/iso_tpm_${cond}.txt"

        REF_PSI_PREFIX="${EVENTS_DIR}/events_${SUPPA_REF_CONDITION}"
        REF_PSI="${REF_PSI_PREFIX}.psi"

        COND_PSI_PREFIX="${EVENTS_DIR}/events_${cond}"
        COND_PSI="${COND_PSI_PREFIX}.psi"
    fi

    require_file "${REF_ISO}"
    require_file "${COND_ISO}"

    if [[ ! -s "${REF_PSI}" ]]; then
        echo "[SUPPA2] Calculating PSI for reference ${SUPPA_REF_CONDITION} vs ${cond}"

        suppa.py psiPerEvent \
            -i "${IOE_MASTER}" \
            -e "${REF_ISO}" \
            -o "${REF_PSI_PREFIX}"

        require_file "${REF_PSI}"
    else
        echo "[SUPPA2] Existing reference PSI found; skipping ${REF_PSI}"
    fi

    if [[ ! -s "${COND_PSI}" ]]; then
        echo "[SUPPA2] Calculating PSI for ${cond} vs ${SUPPA_REF_CONDITION}"

        suppa.py psiPerEvent \
            -i "${IOE_MASTER}" \
            -e "${COND_ISO}" \
            -o "${COND_PSI_PREFIX}"

        require_file "${COND_PSI}"
    else
        echo "[SUPPA2] Existing condition PSI found; skipping ${COND_PSI}"
    fi

    echo "[SUPPA2] Running diffSplice for ${SUPPA_REF_CONDITION} vs ${cond}"

    suppa.py diffSplice \
        -m empirical \
        -gc \
        -s \
        -i "${IOE_MASTER}" \
        -p "${REF_PSI}" "${COND_PSI}" \
        -e "${REF_ISO}" "${COND_ISO}" \
        -o "${OUT_PREFIX}"
done


# -------------------------------
# Mark Complete
# -------------------------------
date -Is >"${DONE_FILE}"

echo "[SUPPA2] Done."
echo "  events:     ${IOE_MASTER}"
echo "  iso_tpm:    ${ISO_TPM_DIR}"
echo "  psi:        ${EVENTS_DIR}"
echo "  diffsplice: ${DIFF_DIR}"
