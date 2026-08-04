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
: "${SALMON_INDEX:?Set SALMON_INDEX in the config}"
: "${STAR_INDEX:?Set STAR_INDEX in the config}"
: "${THREADS:=8}"

metadata="${WORKFLOW_ROOT}/metadata/gse145774_normoxia_0h.tsv"
raw_dir="${OUTPUT_ROOT}/fastq"
salmon_root="${OUTPUT_ROOT}/quant/salmon"
star_root="${OUTPUT_ROOT}/quant/star"
mkdir -p "$salmon_root" "$star_root"

for tool in salmon STAR samtools; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "Required command not found: $tool" >&2
        exit 1
    }
done

[[ -d "$SALMON_INDEX" ]] || { echo "Salmon index not found: $SALMON_INDEX" >&2; exit 1; }
[[ -d "$STAR_INDEX" ]] || { echo "STAR index not found: $STAR_INDEX" >&2; exit 1; }

tail -n +2 "$metadata" |
while IFS=$'\t' read -r sample_id population individual gsm run_accession stress time library_layout read_length; do
    r1="${raw_dir}/${sample_id}_1.fastq.gz"
    r2="${raw_dir}/${sample_id}_2.fastq.gz"
    [[ -s "$r1" && -s "$r2" ]] || {
        echo "FASTQ files missing for $sample_id; run download_fastq.sh first." >&2
        exit 1
    }

    salmon_out="${salmon_root}/${sample_id}"
    if [[ -s "${salmon_out}/quant.sf" ]]; then
        echo "[salmon] present: $sample_id"
    else
        echo "[salmon] $sample_id"
        mkdir -p "$salmon_out"
        salmon quant \
            -i "$SALMON_INDEX" \
            --libType A \
            -1 "$r1" \
            -2 "$r2" \
            --threads "$THREADS" \
            --rangeFactorizationBins 4 \
            --gcBias \
            --seqBias \
            --recoverOrphans \
            --validateMappings \
            -o "$salmon_out"
    fi

    star_out="${star_root}/${sample_id}"
    bam="${star_out}/Aligned.sortedByCoord.out.bam"
    if [[ -s "$bam" ]] && samtools quickcheck "$bam"; then
        echo "[STAR] present: $sample_id"
    else
        echo "[STAR] $sample_id"
        mkdir -p "$star_out"
        STAR \
            --genomeDir "$STAR_INDEX" \
            --readFilesIn "$r1" "$r2" \
            --readFilesCommand zcat \
            --runThreadN "$THREADS" \
            --twopassMode Basic \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS nM \
            --outSAMunmapped Within \
            --outFileNamePrefix "${star_out}/"
        samtools index -@ "$THREADS" "$bam"
        samtools quickcheck "$bam"
    fi
done

echo "[quantify] completed: ${OUTPUT_ROOT}/quant"
