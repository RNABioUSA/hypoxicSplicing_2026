#!/usr/bin/env bash
set -euo pipefail

# Resolve repo root based on this file's location
UTILS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${UTILS_DIR}/../.." && pwd)"

CONFIG_ENV="${PROJECT_ROOT}/config/environment/analysis_environment.env"
[[ -s "${CONFIG_ENV}" ]] || {
    echo "ERROR: Missing ${CONFIG_ENV}" >&2
    exit 1
}

# shellcheck disable=SC1090
source "${CONFIG_ENV}"

# Common repo-standard locations
SAMPLES_TSV="${PROJECT_ROOT}/config/metadata/samples.tsv"
READS_ROOT="${PROJECT_ROOT}/reads"
RESULTS_DIR="${PROJECT_ROOT}/results"

GENOME_DIR="${PROJECT_ROOT}/resources/genome"
GENOME_FA="${GENOME_DIR}/${GENOME_FA}"
GENOME_GTF="${GENOME_DIR}/${GENOME_GTF}"

mkdir -p "${RESULTS_DIR}"
