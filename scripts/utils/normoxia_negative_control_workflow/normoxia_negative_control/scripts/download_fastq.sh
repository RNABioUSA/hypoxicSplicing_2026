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
: "${THREADS:=8}"

metadata="${WORKFLOW_ROOT}/metadata/gse145774_normoxia_0h.tsv"
raw_dir="${OUTPUT_ROOT}/fastq"
tmp_dir="${OUTPUT_ROOT}/sra_tmp"
mkdir -p "$raw_dir" "$tmp_dir"

for tool in fasterq-dump; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "Required command not found: $tool" >&2
        exit 1
    }
done

compressor=gzip
if command -v pigz >/dev/null 2>&1; then
    compressor=pigz
fi

tail -n +2 "$metadata" |
while IFS=$'\t' read -r sample_id population individual gsm run_accession stress time library_layout read_length; do
    r1="${raw_dir}/${sample_id}_1.fastq.gz"
    r2="${raw_dir}/${sample_id}_2.fastq.gz"

    if [[ -s "$r1" && -s "$r2" ]]; then
        echo "[download] present: $sample_id ($run_accession)"
        continue
    fi

    echo "[download] $sample_id ($run_accession)"
    fasterq-dump \
        --split-files \
        --threads "$THREADS" \
        --temp "$tmp_dir" \
        --outdir "$raw_dir" \
        "$run_accession"

    run_r1="${raw_dir}/${run_accession}_1.fastq"
    run_r2="${raw_dir}/${run_accession}_2.fastq"
    [[ -s "$run_r1" && -s "$run_r2" ]] || {
        echo "Expected paired FASTQ files were not produced for $run_accession" >&2
        exit 1
    }

    mv "$run_r1" "${raw_dir}/${sample_id}_1.fastq"
    mv "$run_r2" "${raw_dir}/${sample_id}_2.fastq"

    if [[ "$compressor" == pigz ]]; then
        pigz -p "$THREADS" "${raw_dir}/${sample_id}_1.fastq" "${raw_dir}/${sample_id}_2.fastq"
    else
        gzip "${raw_dir}/${sample_id}_1.fastq" "${raw_dir}/${sample_id}_2.fastq"
    fi

    gzip -t "$r1"
    gzip -t "$r2"
done

echo "[download] completed: $raw_dir"
