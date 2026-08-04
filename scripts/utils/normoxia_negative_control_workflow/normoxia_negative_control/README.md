# Normoxia negative-control analysis

## Recommendation

Do **not** add an external-control label to the main manuscript `samples.tsv` and test the manuscript normoxia samples against an external study as though `condition` were the only model term. In that comparison, study, donor source, culture medium, library preparation, sequencing run, and “condition” are perfectly confounded. Thousands of calls would therefore be unsurprising but would not estimate the background rate of alternative splicing without intervention.

The defensible primary analysis is a separate, internally controlled null analysis:

1. Reprocess untreated HUVEC libraries from a single external study using the same GRCh38.p14/GENCODE v45 annotation and the same Salmon/STAR settings used in the manuscript.
2. Create repeated pseudo-comparisons among untreated samples.
3. Run the manuscript's DESeq2, DEXSeq, DRIMSeq, and SUPPA2 procedures and report the distribution of significant calls, tested denominators, and significant fractions.
4. Use joint PCA of the manuscript and external controls only as sensitivity/QC, with study shown explicitly. Do not interpret the cross-study axis as a hypoxia effect.

The selected dataset is Xin et al. ([Nature Communications 2020](https://doi.org/10.1038/s41467-020-18638-8); [GSE145774](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145774), PRJNA608195/SRP250432). It contains ten independent untreated 0-hour normoxia HUVEC RNA-seq samples: five Tibetan and five Han donors. The libraries are 150-bp paired-end TruSeq mRNA-seq on HiSeq 4000, making them the closest technical match among the studies named by the reviewer.

The five-day normoxia series used in Supplementary Fig. 2B of Xin et al. is not among the 50 RNA-seq libraries listed in the GEO/SRA deposit. This workflow therefore uses the ten deposited 0-hour normoxia libraries. Each pseudo-comparison uses four samples per group, with two Tibetan and two Han donors in each group; one donor from each population is withheld. This balances ancestry without requiring a covariate that SUPPA2 cannot model. The default is 20 unique, reproducible splits.

## Why the other cited studies are weaker choices

| Study | Relevant data | Limitation for this analysis |
|---|---|---|
| [Xin et al. 2020](https://doi.org/10.1038/s41467-020-18638-8) | GSE145774; ten independent 0-hour normoxia HUVEC samples; poly(A), 150-bp paired-end | Best available choice. The paper's additional normoxia time series is not in the listed GEO/SRA RNA-seq samples. |
| [Moreau et al. 2018](https://doi.org/10.3389/fcvm.2018.00159) | GSE118530; two normoxia HUVEC RNA-seq samples | Only two controls; 50-cycle single-end, rRNA-depleted libraries. A 1-vs-1 null test cannot estimate biological dispersion, and the data are less suitable for isoform analysis. |
| [Weigand et al. 2012](https://doi.org/10.1371/journal.pone.0042697) | HUVEC normoxia/hypoxia profiled with Affymetrix Human Exon 1.0 ST arrays | Not RNA-seq; cannot be reprocessed with Salmon, DEXSeq, DRIMSeq, or SUPPA2. |
| [Niskanen et al. 2018](https://doi.org/10.1093/nar/gkx1214) | GSE94872 plus referenced HUVEC expression data | The hypoxia transcription assay is GRO-seq and the main study focuses on TCC/Hi-C; deposited reads are short single-end. This is not comparable steady-state, paired-end mRNA-seq for isoform usage. |

## Relationship to the uploaded manuscript scripts

This workflow mirrors the manuscript analysis choices where they affect significance:

- DESeq2: gene-level BH FDR below 0.10.
- DEXSeq: paired-end, unstranded `summarizeOverlaps`, multi-overlap and multimapping reads retained, exon-bin BH FDR and `perGeneQValue` at 0.10.
- DRIMSeq: `dtuScaledTPM` counts, the manuscript's H1/C1 filters (`3, 6, 10, 10`), and stageR OFDR 0.10.
- SUPPA2: empirical `diffSplice`, `-gc`, and `-s`.

The uploaded `load_suppa.R` treats an unadjusted SUPPA2 p-value below 0.10 as significant. That is not directly comparable to the 0.10 FDR/OFDR thresholds used by the other tools. The companion analysis therefore reports both:

- `SUPPA2_rawP`, reproducing the current manuscript rule; and
- `SUPPA2_BH`, applying BH correction across tested events, recommended for the primary cross-tool summary.

DEXSeq needs genome-aligned BAM files. Adjusting only the Salmon metadata is therefore insufficient for the requested four-tool analysis. The external FASTQs are quantified with Salmon and aligned with STAR once; all pseudo-comparisons reuse those results.

## Files

- `metadata/gse145774_normoxia_0h.tsv`: verified GEO/SRA sample-to-run mapping.
- `config/negative_control.env.example`: path and run settings.
- `scripts/download_fastq.sh`: downloads and compresses the ten paired FASTQ libraries.
- `scripts/quantify_align.sh`: runs Salmon using the manuscript options and creates STAR BAMs for DEXSeq.
- `scripts/run_null_statistics.R`: creates balanced pseudo-groups and runs DESeq2, DRIMSeq, and DEXSeq.
- `scripts/run_suppa_null.sh`: runs SUPPA2 for each split.
- `scripts/summarize_null.R`: adds BH-adjusted SUPPA2 results, tested denominators, gene-level intersections, summaries, and a count-distribution plot.
- `scripts/run_all.sh`: executes the full workflow.

## Dependencies

Command-line tools:

```text
SRA Toolkit (fasterq-dump)
Salmon
STAR
samtools
SUPPA2 (suppa.py)
pigz (optional; gzip is used otherwise)
```

R/Bioconductor packages:

```r
BiocManager::install(c(
  "tximport", "DESeq2", "DRIMSeq", "stageR", "DEXSeq",
  "rtracklayer", "GenomicFeatures", "GenomicAlignments",
  "Rsamtools", "SummarizedExperiment", "S4Vectors", "BiocParallel"
))
```

## Run

Copy the workflow into the manuscript project or another analysis location, then:

```bash
cp config/negative_control.env.example config/negative_control.env
# Edit all paths in config/negative_control.env.

bash scripts/run_all.sh config/negative_control.env
```

For a short pilot, set `N_PERM=2`; for the revision, use at least 20. These ten libraries contain roughly 35–44 million read pairs each. Plan for substantial download, FASTQ, BAM, and temporary storage (approximately 300–500 GB of scratch space is prudent).

## Main outputs

Under `OUTPUT_ROOT/analysis`:

- `null_gene_count_summary.tsv`: median, IQR, and range of DEGs/AS genes and the three-tool splicing intersection.
- `null_gene_counts_by_permutation.tsv`: one row per pseudo-comparison.
- `all_tool_null_summary.tsv`: gene-, transcript-, exon-bin-, and event-level counts, tested denominators, and fractions.
- `all_tool_results_by_permutation.tsv`: long-format results for every method and split.
- `all_significant_gene_ids.tsv`: significant genes for overlap checks.
- `null_gene_count_distributions.png`: concise supplementary figure candidate.
- `permutations/perm_*/`: full DESeq2, DRIMSeq, and DEXSeq result tables.
- `suppa/perm_*/`: SUPPA2 dPSI, p-value, and BH-adjusted event tables.

The manuscript should report distributions, not a single favorable split. A compact result sentence is:

> Across 20 population-balanced pseudo-comparisons of untreated HUVEC samples, the median [IQR; range] numbers of significant genes were **[DESeq2]**, **[DEXSeq]**, **[DRIMSeq]**, and **[SUPPA2-BH]**, while the median three-way intersection of alternatively spliced genes was **[intersection]**. The corresponding hypoxia contrasts yielded **[insert manuscript counts]**.

## Draft response to the reviewer

> We thank the reviewer for raising the important question of the background number of differential-expression and alternative-splicing calls expected without an intervention. We agree that a direct comparison between our normoxia samples and public controls would be confounded by study-specific culture, donor, library-preparation, and sequencing effects. We therefore performed a separate, internally controlled null analysis using the ten independent untreated (0 h, normoxia) HUVEC RNA-seq samples from Xin et al. (GSE145774), the closest of the suggested datasets in cell type and sequencing design. All reads were reprocessed with the same GRCh38.p14/GENCODE v45 reference and Salmon/STAR workflow used for our data. We generated 20 reproducible population-balanced 4-versus-4 pseudo-comparisons and applied DESeq2, DEXSeq, DRIMSeq/stageR, and SUPPA2 using the manuscript thresholds. Across these comparisons, the median [IQR; range] numbers of significant genes were **[DESeq2 result]**, **[DEXSeq result]**, **[DRIMSeq result]**, and **[SUPPA2-BH result]**; the median three-method intersection for alternative splicing was **[intersection result]**. These background counts were substantially lower than those observed after 1, 3, and 24 h hypoxia (**[insert corresponding observed counts/fractions]**). We have added the analysis as Supplementary Fig. **[X]** and Supplementary Table **[Y]**, report both the number tested and fraction significant for each method, and clarified the limitations of inferring pericellular oxygen kinetics without direct dissolved-oxygen measurements.

If the null and hypoxia distributions overlap materially, replace “substantially lower” with a neutral description and temper the manuscript's claim. This control assesses background analysis calls; it does not itself prove that the medium reached the intended oxygen concentration within one hour.
