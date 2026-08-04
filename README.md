# Hypoxic Splicing in Human Endothelial Cells

This repository contains the computational analysis code accompanying a manuscript examining transcriptional and alternative-splicing responses to hypoxia in human umbilical vein endothelial cells (HUVECs).

The analyses include:

- transcript quantification with Salmon
- differential gene-expression analysis with DESeq2
- differential exon-usage analysis with DEXSeq
- differential transcript-usage analysis with DRIMSeq
- event-level alternative-splicing analysis with SUPPA2
- functional enrichment analyses
- generation of manuscript figures and supplemental data tables

## Repository organization

- `scripts/analysis/` — primary analysis workflows
- `scripts/plots/` — manuscript figure-generation code
- `scripts/prep/` — preprocessing and feature-counting utilities
- `scripts/utils/` — shared helper functions and supplemental table generation
- `config/` — analysis configuration and sample metadata
- `resources/` — reference resources used by the workflows
- `results/` — locally generated analysis outputs

## Analysis workflow

The accompanying R Notebook serves as the primary guide to the analysis. It documents the execution order for the command-line and R workflows and provides the code used to generate manuscript figures and supplemental tables.

## Data and reference files

Raw sequencing data, genome reference files, and other large inputs are not distributed through this repository. These resources must be obtained separately and placed in the locations expected by the analysis configuration and sample metadata.

## Repository version

This repository contains the code provided for journal review. A manuscript citation, data accession, and archived release identifier will be added when available.

## License

This repository is distributed under the MIT License. See `LICENSE` for details.
