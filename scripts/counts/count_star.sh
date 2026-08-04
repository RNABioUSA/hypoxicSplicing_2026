#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Script: count_star.sh
#
# Purpose
# -------
# Perform splice-aware genome alignment with STAR for all samples listed in the
# project sample sheet, producing coordinate-sorted BAM files suitable for
# downstream exon-level counting for DEXSeq.
#
# The script:
#
#   1. Builds a STAR genome index (if it does not already exist).
#   2. Iterates through the unified sample metadata table.
#   3. Aligns each sample with STAR using ENCODE-like RNA-seq parameters.
#   4. Indexes the resulting BAM files for downstream counting / QC.
#
# The script avoids recomputing work when outputs already exist.
#
# Workflow Overview
# -----------------
# FASTQ reads
#      │
#      ▼
# STAR genomeGenerate
#      │
#      ▼
# STAR alignReads
#      │
#      ▼
# Aligned.sortedByCoord.out.bam
#      │
#      ▼
# samtools index
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
# resources/star_index/<GENOME_BUILD>/
#   STAR genome index
#
# results/counts/star/
#   <condition>/<sample>/
#       <sample>.Aligned.sortedByCoord.out.bam
#       <sample>.Aligned.sortedByCoord.out.bam.bai
#       <sample>.Log.final.out
#       <sample>.SJ.out.tab
#   <condition>/
#       <condition>.Aligned.sortedByCoord.out.bam
#       <condition>.Aligned.sortedByCoord.out.bam.bai
#       <condition>.SJ.raw.out.bed
#       <condition>.SJ.cpm.out.bed
#       <condition>.bw
#
# Dependencies
# ------------
# External tools (expected in PATH via conda environment):
#
#   STAR
#   samtools
#   deeptools
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
#       "$PROJECT_ROOT/scripts/counts/count_star.sh"
#
# Notes for DEXSeq
# ----------------
# - These alignments are intended for downstream exon/bin counting.
# - Strand handling for DEXSeq is controlled during counting
#   (e.g. dexseq_count.py / summarizeOverlaps), not by STAR itself.
# - The optional XS-tag behavior (--outSAMstrandField intronMotif) is not
#   required for DEXSeq and is therefore left off by default.
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

star_index_exists() {
    local f
    for f in "${STAR_INDEX_FILES[@]}"; do
        [[ -s "${f}" ]] || return 1
    done
    return 0
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

require_tool STAR
require_tool samtools
require_tool bamCoverage

STAR_INDEX="${PROJECT_ROOT}/resources/star_index"
OUTDIR="${RESULTS_DIR}/counts/star"

mkdir -p "${OUTDIR}"

DONE_FILE="${OUTDIR}/.star_align.done"

# -------------------------------
# Set Paths Paramters
# -------------------------------

# STAR recommends sjdbOverhang = readLength - 1.
SJDB_OVERHANG=149

# -------------------------------
# Merged Paramters
# -------------------------------
BIGWIG_BINSIZE=10
BIGWIG_NORM="CPM"

# -------------------------------
# Global Skip Check
# -------------------------------
if [[ -s "${DONE_FILE}" ]]; then
    echo "[STAR] Alignment already complete (${DONE_FILE}); skipping."
    exit 0
fi

# -------------------------------
# Build STAR Index (skip if exists)
# -------------------------------
STAR_INDEX_FILES=(
    "${STAR_INDEX}/Genome"
    "${STAR_INDEX}/SA"
    "${STAR_INDEX}/SAindex"
    "${STAR_INDEX}/chrLength.txt"
    "${STAR_INDEX}/chrName.txt"
    "${STAR_INDEX}/chrStart.txt"
    "${STAR_INDEX}/sjdbInfo.txt"
    "${STAR_INDEX}/genomeParameters.txt"
)

if star_index_exists; then
    echo "[STAR] Found existing index: ${STAR_INDEX}"
else
    echo "[STAR] STAR index not found; building ${GENOME_BUILD} index..."

    require_file "${GENOME_FA}"
    require_file "${GENOME_GTF}"

    mkdir -p "${STAR_INDEX}"

    STAR \
        --runThreadN "${THREADS}" \
        --runMode genomeGenerate \
        --genomeDir "${STAR_INDEX}" \
        --genomeFastaFiles "${GENOME_FA}" \
        --sjdbGTFfile "${GENOME_GTF}" \
        --sjdbOverhang "${SJDB_OVERHANG}"

    echo "[STAR] Index build complete."
fi

# -------------------------------
# Set STAR Parameters
# -------------------------------
# ENCODE-like STAR settings as described in STAR documentation.
STAR_COMMON=(
    --genomeDir "${STAR_INDEX}"
    --readFilesCommand zcat
    --outSAMtype BAM SortedByCoordinate
    --outSAMunmapped Within
    --outFilterType BySJout
    --outSAMattributes NH HI AS NM MD
    --outFilterMultimapNmax 20
    --outFilterMismatchNmax 999
    --outFilterMismatchNoverReadLmax 0.04
    --alignIntronMin 20
    --alignIntronMax 1000000
    --alignMatesGapMax 1000000
    --alignSJoverhangMin 8
    --alignSJDBoverhangMin 1
    --sjdbScore 1
    --twopassMode Basic
)

# -------------------------------
# Align per sample from samples.tsv
# Required columns:
#   sample_id, condition, read1_rel, read2_rel, layout, stranded, exclude
# -------------------------------
tail -n +2 "${SAMPLES_TSV}" |
    while IFS=$'\t' read -r sample cond rep r1_rel r2_rel layout stranded excl; do

        [[ -z "${sample}" ]] && continue
        excl="${excl//$'\r'/}"

        if [[ "${excl}" == "TRUE" ]]; then
            echo "[star] Skipping excluded sample ${sample}"
            continue
        fi

        R1="${READS_ROOT}/${r1_rel}"
        R2="${READS_ROOT}/${r2_rel}"

        require_file "${R1}"
        require_file "${R2}"

        sample_out="${OUTDIR}/${cond}/${sample}"
        prefix="${sample_out}/${sample}."

        bam="${prefix}Aligned.sortedByCoord.out.bam"
        bai="${bam}.bai"
        final_log="${prefix}Log.final.out"

        mkdir -p "${sample_out}"

        if [[ -s "${bam}" && -s "${bai}" && -s "${final_log}" ]]; then
            echo "[star] Existing alignment found; skipping ${sample}"
            continue
        fi

        echo "[star] Aligning ${sample} (${cond})"

        STAR \
            --runThreadN "${THREADS}" \
            "${STAR_COMMON[@]}" \
            --readFilesIn "${R1}" "${R2}" \
            --outFileNamePrefix "${prefix}" \
            --outSAMattrRGline "ID:${sample}" "SM:${sample}"

        samtools index -@ "${THREADS}" "${bam}"
    done
    
# -------------------------------
# Pool replicate BAMs by condition for visualization
# Generate:
#   - pooled BAM
#   - pooled BAM index
#   - pooled raw STAR junction table
#   - pooled CPM-normalized STAR junction table
#   - pooled CPM-normalized bigWig
# -------------------------------
echo "[STAR] Pooling replicate BAMs by condition and generating sashimi files"

mapfile -t CONDITIONS < <(
    tail -n +2 "${SAMPLES_TSV}" |
    awk -F '\t' '$8 != "TRUE" && $2 != "" {print $2}' |
    sort -u
)

for cond in "${CONDITIONS[@]}"; do
    cond_dir="${OUTDIR}/${cond}"
    mkdir -p "${cond_dir}"

    pooled_bam="${cond_dir}/${cond}.Aligned.sortedByCoord.out.bam"
    pooled_bai="${pooled_bam}.bai"
    pooled_bw="${cond_dir}/${cond}.bw"
    pooled_sj_raw="${cond_dir}/${cond}.SJ.raw.out.tab"
    pooled_sj_cpm="${cond_dir}/${cond}.SJ.cpm.out.tab"

    mapfile -t bam_list < <(
        find "${cond_dir}" -mindepth 2 -maxdepth 2 -type f -name "*.Aligned.sortedByCoord.out.bam" | sort
    )

    mapfile -t sj_list < <(
        find "${cond_dir}" -mindepth 2 -maxdepth 2 -type f -name "*.SJ.out.tab" | sort
    )

    if [[ "${#bam_list[@]}" -eq 0 ]]; then
        echo "[STAR] No replicate BAMs found for condition ${cond}; skipping pooled outputs"
        continue
    fi

    if [[ "${#sj_list[@]}" -eq 0 ]]; then
        echo "[STAR] No replicate SJ.out.tab files found for condition ${cond}; skipping pooled junction outputs"
        continue
    fi

    echo "[STAR] Condition ${cond}: found ${#bam_list[@]} replicate BAM(s)"
    echo "[STAR] Condition ${cond}: found ${#sj_list[@]} replicate SJ.out.tab file(s)"

    if [[ ! -s "${pooled_bam}" ]]; then
        echo "[STAR] Merging replicate BAMs for ${cond}"
        samtools merge -@ "${THREADS}" -f "${pooled_bam}" "${bam_list[@]}"
    else
        echo "[STAR] Existing pooled BAM found for ${cond}; skipping merge"
    fi

    if [[ ! -s "${pooled_bai}" ]]; then
        echo "[STAR] Indexing pooled BAM for ${cond}"
        samtools index -@ "${THREADS}" "${pooled_bam}"
    else
        echo "[STAR] Existing pooled BAM index found for ${cond}; skipping"
    fi

    if [[ ! -s "${pooled_sj_raw}" ]]; then
        echo "[STAR] Merging replicate SJ.out.tab files for ${cond}"

        awk 'BEGIN{OFS="\t"}
        {
            key = $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $9
            uniq[key] += $7
            multi[key] += $8
        }
        END {
            for (k in uniq) {
                split(k, a, OFS)
                print a[1], a[2], a[3], a[4], a[5], a[6], uniq[k], multi[k], a[7]
            }
        }' "${sj_list[@]}" | sort -k1,1 -k2,2n -k3,3n > "${pooled_sj_raw}"
    else
        echo "[STAR] Existing pooled raw SJ.out.tab found for ${cond}; skipping merge"
    fi

    if [[ ! -s "${pooled_sj_cpm}" ]]; then
        echo "[STAR] Normalizing pooled SJ.out.tab to CPM for ${cond}"

        mapped_reads=$(samtools view -c -F 4 "${pooled_bam}")
        [[ "${mapped_reads}" -gt 0 ]] || die "Mapped read count is zero for ${pooled_bam}"

        awk -v OFS="\t" -v M="${mapped_reads}" '
        {
            $7 = ($7 / M) * 1000000
            $8 = ($8 / M) * 1000000
            print
        }' "${pooled_sj_raw}" > "${pooled_sj_cpm}"
    else
        echo "[STAR] Existing pooled CPM-normalized SJ.out.tab found for ${cond}; skipping"
    fi

    if [[ ! -s "${pooled_bw}" ]]; then
        echo "[STAR] Generating pooled bigWig for ${cond}"
        bamCoverage \
            -b "${pooled_bam}" \
            -o "${pooled_bw}" \
            --outFileFormat bigwig \
            --binSize "${BIGWIG_BINSIZE}" \
            --normalizeUsing "${BIGWIG_NORM}" \
            --exactScaling \
            --numberOfProcessors "${THREADS}"
    else
        echo "[STAR] Existing pooled bigWig found for ${cond}; skipping"
    fi
done

date -Is >"${DONE_FILE}"

echo "[STAR] Done."
echo "  index:      ${STAR_INDEX}"
echo "  alignments: ${OUTDIR}/<condition>/<sample>/<sample>.Aligned.sortedByCoord.out.bam"
