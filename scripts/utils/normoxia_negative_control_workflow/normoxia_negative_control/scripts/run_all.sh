#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 CONFIG_ENV" >&2
    exit 2
}

[[ $# -eq 1 ]] || usage
config_file=$1
[[ -s "$config_file" ]] || { echo "Config not found: $config_file" >&2; exit 1; }

# shellcheck disable=SC1090
source "$config_file"

: "${WORKFLOW_ROOT:?Set WORKFLOW_ROOT in the config}"
: "${OUTPUT_ROOT:?Set OUTPUT_ROOT in the config}"
: "${GTF:?Set GTF in the config}"
: "${N_PERM:=20}"
: "${SEED:=20260727}"
: "${ALPHA:=0.10}"
: "${THREADS:=8}"

"${WORKFLOW_ROOT}/scripts/download_fastq.sh" "$config_file"
"${WORKFLOW_ROOT}/scripts/quantify_align.sh" "$config_file"

Rscript "${WORKFLOW_ROOT}/scripts/run_null_statistics.R" \
    --metadata "${WORKFLOW_ROOT}/metadata/gse145774_normoxia_0h.tsv" \
    --quant-root "${OUTPUT_ROOT}/quant" \
    --gtf "$GTF" \
    --out "${OUTPUT_ROOT}/analysis" \
    --n-perm "$N_PERM" \
    --seed "$SEED" \
    --alpha "$ALPHA" \
    --workers "$THREADS"

"${WORKFLOW_ROOT}/scripts/run_suppa_null.sh" "$config_file"

Rscript "${WORKFLOW_ROOT}/scripts/summarize_null.R" \
    --analysis-root "${OUTPUT_ROOT}/analysis" \
    --alpha "$ALPHA"

echo "[workflow] completed: ${OUTPUT_ROOT}/analysis"
