#!/usr/bin/env bash
set -euo pipefail

echo "🔎 Checking environment..."

fail=0

need_cmd() {
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "❌ Missing: $1"
        fail=1
    else
        echo "✅ Found: $1 ($(command -v "$1"))"
    fi
}

echo
echo "Checking core bioinformatics tools:"
need_cmd salmon
need_cmd STAR
need_cmd samtools
need_cmd featureCounts
need_cmd gffread
need_cmd infer_experiment.py

echo
echo "Checking compression tools:"
need_cmd gzip
need_cmd gunzip
need_cmd pigz

echo
echo "Checking R availability:"
need_cmd R
need_cmd Rscript
need_cmd Rscript

echo
echo "Checking R packages (via Rscript)..."

Rscript - <<'EOF' || fail=1
required <- c(
  "tximeta",
  "DESeq2",
  "DEXSeq",
  "GenomicAlignments",
  "Rsamtools",
  "rtracklayer",
  "SummarizedExperiment",
  "dplyr",
  "ggplot2",
  "readr",
  "stringr"
)

missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]

if (length(missing)) {
  cat("❌ Missing R packages:\n")
  print(missing)
  quit(status = 1)
} else {
  cat("✅ All required R packages found\n")
}
EOF

echo
if [[ "$fail" -eq 1 ]]; then
    echo "❌ Environment check FAILED"
    exit 1
else
    echo "✅ Environment check PASSED"
fi
