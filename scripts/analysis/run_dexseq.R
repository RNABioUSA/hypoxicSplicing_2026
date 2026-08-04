PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
if (!nzchar(PROJECT_ROOT)) stop("PROJECT_ROOT env var not set.")

HELPERS_FILE <- file.path(PROJECT_ROOT, "scripts/utils/helpers.R")
if (!file.exists(HELPERS_FILE)) stop("helpers.R not found at: ", HELPERS_FILE)
source(HELPERS_FILE)

# -------------------------
# Internal Helpers
# -------------------------

.standardize_dexseq_gene_res <- function(res_df,
                                         timepoint,
                                         contrast,
                                         tool = "dexseq",
                                         annot_df = NULL) {
  out <- as.data.frame(res_df, stringsAsFactors = FALSE)

  if (!("groupID" %in% colnames(out))) {
    stop("DEXSeq gene-level results must contain column 'groupID'.")
  }

  out$gene_id_full <- out$groupID
  out$gene_id <- .strip_ens_version(out$groupID)
  out$ensgene <- out$gene_id

  out$timepoint <- timepoint
  out$contrast <- contrast
  out$comparison <- contrast
  out$tool <- tool

  if (!("padj" %in% colnames(out))) {
    out$padj <- NA_real_
  }

  out$neglog10_padj <- ifelse(
    is.na(out$padj) | out$padj <= 0,
    NA_real_,
    -log10(out$padj)
  )

  out <- .attach_gene_annot(out, annot_df = annot_df)
  out
}

.standardize_dexseq_exon_res <- function(res_df,
                                         timepoint,
                                         contrast,
                                         tool = "dexseq",
                                         annot_df = NULL) {
  out <- as.data.frame(res_df, stringsAsFactors = FALSE)

  if (!all(c("groupID", "featureID") %in% colnames(out))) {
    stop("DEXSeq exon-level results must contain columns 'groupID' and 'featureID'.")
  }

  out$gene_id_full <- out$groupID
  out$gene_id <- .strip_ens_version(out$groupID)
  out$ensgene <- out$gene_id

  out$exon_id <- out$featureID
  out$feature_id <- paste(out$groupID, out$featureID, sep = ":")

  out$timepoint <- timepoint
  out$contrast <- contrast
  out$comparison <- contrast
  out$tool <- tool

  if (!("padj" %in% colnames(out))) {
    out$padj <- NA_real_
  }

  out$neglog10_padj <- ifelse(
    is.na(out$padj) | out$padj <= 0,
    NA_real_,
    -log10(out$padj)
  )

  out <- .attach_gene_annot(out, annot_df = annot_df)
  out
}

.make_dexseq_summary <- function(gene_sig_list,
                                 exon_sig_list,
                                 padj_cutoff = 0.10) {
  tps <- names(gene_sig_list)
  do.call(rbind, lapply(tps, function(tp) {
    data.frame(
      timepoint = tp,
      padj_cutoff = padj_cutoff,
      n_sig_genes = nrow(gene_sig_list[[tp]]),
      n_sig_exons = nrow(exon_sig_list[[tp]]),
      stringsAsFactors = FALSE
    )
  }))
}

.write_dexseq_xlsx <- function(gene_sig_list,
                               exon_sig_list,
                               summary_tbl,
                               out_xlsx) {
  .check_pkg("openxlsx")

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Summary")
  openxlsx::writeData(wb, "Summary", summary_tbl)

  for (tp in names(gene_sig_list)) {
    gene_sheet <- paste0(tp, "_genes")
    exon_sheet <- paste0(tp, "_exons")

    openxlsx::addWorksheet(wb, gene_sheet)
    if (nrow(gene_sig_list[[tp]]) == 0) {
      openxlsx::writeData(
        wb, gene_sheet,
        data.frame(note = "No significant genes at this cutoff.")
      )
    } else {
      openxlsx::writeData(wb, gene_sheet, gene_sig_list[[tp]])
    }

    openxlsx::addWorksheet(wb, exon_sheet)
    if (nrow(exon_sig_list[[tp]]) == 0) {
      openxlsx::writeData(
        wb, exon_sheet,
        data.frame(note = "No significant exons at this cutoff.")
      )
    } else {
      openxlsx::writeData(wb, exon_sheet, exon_sig_list[[tp]])
    }
  }

  dir.create(dirname(out_xlsx), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
  invisible(TRUE)
}

.read_dexseq_bam_coldata <- function(samples_tsv,
                                     results_dir = .get_results_dir(),
                                     contrast_var = "condition",
                                     sample_col = "sample_id") {
  x <- utils::read.delim(
    samples_tsv,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required <- c("sample_id", "condition", "exclude")
  missing_cols <- setdiff(required, colnames(x))
  if (length(missing_cols) > 0) {
    stop("samples.tsv missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  x <- x[is.na(x$exclude) | x$exclude != "TRUE", , drop = FALSE]

  x$names <- x[[sample_col]]
  x$bam <- file.path(
    results_dir,
    "counts",
    "star",
    x[[contrast_var]],
    x$names,
    paste0(x$names, ".Aligned.sortedByCoord.out.bam")
  )

  missing_bams <- x$bam[!file.exists(x$bam)]
  if (length(missing_bams) > 0) {
    stop(
      "Missing BAM files. Run STAR first.\n",
      paste(missing_bams, collapse = "\n")
    )
  }

  x
}

.read_dexseq_htseq_coldata <- function(samples_tsv,
                                       results_dir = .get_results_dir(),
                                       contrast_var = "condition",
                                       sample_col = "sample_id") {
  x <- utils::read.delim(
    samples_tsv,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required <- c("sample_id", "condition", "exclude")
  missing_cols <- setdiff(required, colnames(x))
  if (length(missing_cols) > 0) {
    stop("samples.tsv missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  x <- x[is.na(x$exclude) | x$exclude != "TRUE", , drop = FALSE]

  x$names <- x[[sample_col]]
  x$count_file <- file.path(
    results_dir,
    "counts",
    "htseq",
    paste0(x$names, ".dexseq.txt")
  )

  missing_counts <- x$count_file[!file.exists(x$count_file)]
  if (length(missing_counts) > 0) {
    stop(
      "Missing DEXSeq HTSeq count files. Run count_htseq.sh first.\n",
      paste(missing_counts, collapse = "\n")
    )
  }

  x
}

.count_dexseq_summarizeOverlaps <- function(bam_files,
                                            flattened_annotation,
                                            sample_table,
                                            paired_end = TRUE,
                                            ignore_strand = TRUE,
                                            fragments = TRUE,
                                            count_multi_exon_overlaps = TRUE,
                                            include_multimappers = TRUE,
                                            yield_size = 1000000,
                                            BPPARAM = BiocParallel::SerialParam()) {
  .check_pkg("Rsamtools")
  .check_pkg("GenomicAlignments")
  .check_pkg("SummarizedExperiment")
  .check_pkg("S4Vectors")
  .check_pkg("BiocParallel")

  bam_list <- Rsamtools::BamFileList(
    lapply(bam_files, function(x) Rsamtools::BamFile(x, yieldSize = yield_size))
  )

  old_param <- BiocParallel::registered()[[1]]
  BiocParallel::register(BPPARAM, default = TRUE)
  on.exit(BiocParallel::register(old_param, default = TRUE), add = TRUE)

  param <- Rsamtools::ScanBamParam(tag = "NH")

  preprocess_reads <- function(reads) {
    if (include_multimappers) {
      return(reads)
    }

    if (inherits(reads, "GAlignmentsList")) {
      u <- unlist(reads, use.names = FALSE)
      nh <- S4Vectors::mcols(u)$NH

      if (is.null(nh)) {
        warning("NH tag not present in unlisted GAlignmentsList; multimapper filtering was not applied.")
        return(reads)
      }

      nh_by_fragment <- relist(nh, reads)

      keep <- vapply(
        nh_by_fragment,
        function(z) length(z) > 0 && all(!is.na(z) & z == 1),
        logical(1)
      )

      return(reads[keep])
    }

    nh <- S4Vectors::mcols(reads)$NH

    if (is.null(nh)) {
      warning("NH tag not present in reads; multimapper filtering was not applied.")
      return(reads)
    }

    reads[!is.na(nh) & nh == 1]
  }

  se <- GenomicAlignments::summarizeOverlaps(
    features = flattened_annotation,
    reads = bam_list,
    singleEnd = !paired_end,
    fragments = fragments,
    ignore.strand = ignore_strand,
    inter.feature = !count_multi_exon_overlaps,
    param = param,
    preprocess.reads = preprocess_reads
  )

  SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(sample_table)
  se
}

.count_dexseq_featureCounts <- function(bam_files,
                                        fc_gtf,
                                        paired_end = TRUE,
                                        ignore_strand = TRUE,
                                        count_multi_exon_overlaps = TRUE,
                                        include_multimappers = FALSE,
                                        threads = 1) {
  .check_pkg("Rsubread")

  fc_gtf <- .get_featurecounts_exon_gtf(fc_gtf)

  strand_mode <- if (ignore_strand) 0 else 2

  zz <- file(tempfile(), open = "wt")
  sink(zz, type = "message")

  fc <- Rsubread::featureCounts(
    files = bam_files,
    annot.ext = fc_gtf,
    isGTFAnnotationFile = TRUE,
    GTF.featureType = "exon",
    GTF.attrType = "gene_id",
    useMetaFeatures = FALSE,
    allowMultiOverlap = count_multi_exon_overlaps,
    countMultiMappingReads = include_multimappers,
    isPairedEnd = paired_end,
    countReadPairs = paired_end,
    strandSpecific = strand_mode,
    nthreads = threads
  )

  sink(type = "message")
  close(zz)

  fc
}

.count_dexseq_HTSeq <- function(count_files,
                                sample_table,
                                flattenedfile) {
  .check_pkg("DEXSeq")

  if (is.null(flattenedfile) || !file.exists(flattenedfile)) {
    stop("A valid flattened HTSeq GFF must be supplied.")
  }

  if (length(count_files) != nrow(sample_table)) {
    stop("count_files length must match number of rows in sample_table.")
  }

  clean_files <- vapply(
    count_files,
    .clean_dexseq_htseq_countfile,
    character(1)
  )

  names(clean_files) <- rownames(sample_table)

  DEXSeq::DEXSeqDataSetFromHTSeq(
    countfiles = unname(clean_files),
    sampleData = sample_table,
    design = ~ sample + exon + condition:exon,
    flattenedfile = flattenedfile
  )
}

.clean_dexseq_htseq_countfile <- function(count_file, overwrite = TRUE) {
  clean_file <- sub("\\.txt$", ".clean.txt", count_file)

  if (!file.exists(clean_file) || overwrite) {
    x <- readLines(count_file, warn = FALSE)
    x <- trimws(x)
    x <- x[nzchar(x)]
    x <- x[!grepl("^_", x)]

    # Fix malformed DEXSeq HTSeq IDs:
    x <- sub('^"([^"]+)":"([^"]+)"\t([0-9]+)$', "\\1:\\2\t\\3", x)

    writeLines(x, clean_file)
  }

  clean_file
}

.get_dexseq_annotation_path <- function(resources_dir,
                                        gtf_base,
                                        method,
                                        ext) {
  annot_dir <- file.path(resources_dir, "dexseq_annotation")
  dir.create(annot_dir, recursive = TRUE, showWarnings = FALSE)

  file.path(
    annot_dir,
    paste0(gtf_base, ".flat.", method, ".", ext)
  )
}

.get_dexseq_flattened_annotation <- function(env_file,
                                             resources_dir,
                                             overwrite = FALSE) {
  .check_pkg("GenomicFeatures")

  txdb_info <- .get_txdb(
    env_file = env_file,
    resources_dir = resources_dir,
    overwrite = overwrite
  )

  txdb <- txdb_info$txdb

  flat_rds <- .get_dexseq_annotation_path(
    resources_dir = resources_dir,
    gtf_base = txdb_info$gtf_base,
    method = "summarizeoverlaps",
    ext = "rds"
  )

  if (file.exists(flat_rds) && !overwrite) {
    message("Loading cached DEXSeq flattened annotation: ", flat_rds)
    return(readRDS(flat_rds))
  }

  message("Building DEXSeq flattened annotation from TxDb...")
  flattened_annotation <- GenomicFeatures::exonicParts(
    txdb,
    linked.to.single.gene.only = TRUE
  )

  names(flattened_annotation) <- sprintf(
    "%s:E%03d",
    flattened_annotation$gene_id,
    flattened_annotation$exonic_part
  )

  saveRDS(flattened_annotation, flat_rds)
  message("Cached flattened annotation saved to: ", flat_rds)

  flattened_annotation
}

.get_dexseq_featurecounts_gtf <- function(env_file,
                                          resources_dir,
                                          overwrite = FALSE) {
  ref <- .resolve_ref_from_env(env_file, resources_dir)

  fc_gtf <- .get_dexseq_annotation_path(
    resources_dir = resources_dir,
    gtf_base = ref$gtf_base,
    method = "featurecounts",
    ext = "gtf"
  )

  if (file.exists(fc_gtf) && !overwrite) {
    message("Loading cached DEXSeq featureCounts GTF: ", fc_gtf)
    return(fc_gtf)
  }

  stop("DEXSeq featureCounts GTF not found yet; generate it first.")
}

.get_featurecounts_exon_gtf <- function(fc_gtf) {
  exon_gtf <- sub("\\.gtf$", ".exons.gtf", fc_gtf)

  if (!file.exists(exon_gtf)) {
    gtf <- read.delim(fc_gtf, header = FALSE, stringsAsFactors = FALSE)

    exon_gtf_df <- gtf[gtf$V3 == "exon", ]

    write.table(
      exon_gtf_df,
      exon_gtf,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
  }

  exon_gtf
}

.get_dexseq_htseq_gff <- function(env_file,
                                  resources_dir,
                                  overwrite = FALSE) {
  ref <- .resolve_ref_from_env(env_file, resources_dir)

  gff_file <- .get_dexseq_annotation_path(
    resources_dir = resources_dir,
    gtf_base = ref$gtf_base,
    method = "htseq",
    ext = "gff"
  )

  if (!file.exists(gff_file)) {
    stop(
      "DEXSeq HTSeq flattened GFF not found: ", gff_file, "\n",
      "Run the HTSeq preparation/counting script first."
    )
  }

  gff_file
}

.ensure_dexseq_factors <- function(df, ref_level) {
  df$sample <- factor(df$sample)
  df$condition <- stats::relevel(factor(df$condition), ref = ref_level)
  if ("libType" %in% colnames(df)) {
    df$libType <- factor(df$libType)
  }
  df
}

.get_dexseq_mode_tag <- function(counting_method,
                                 count_multi_exon_overlaps = TRUE,
                                 include_multimappers = FALSE) {
  if (counting_method == "HTSeq") {
    return("HTSeq.multiOverlap.unique")
  }

  if (counting_method == "summarizeOverlaps") {
    return(paste0(
      "summarizeOverlaps",
      if (count_multi_exon_overlaps) ".multiOverlap" else ".strictOverlap",
      if (include_multimappers) ".all" else ".unique"
    ))
  }

  if (counting_method == "featureCounts") {
    return(paste0(
      "featureCounts",
      if (count_multi_exon_overlaps) ".multiOverlap" else ".strictOverlap",
      if (include_multimappers) ".all" else ".unique"
    ))
  }

  stop("Unsupported counting_method: ", counting_method)
}


.DEXSeqDataSetFromFeatureCounts <- function(fc,
                                            sampleData,
                                            design = ~ sample + exon + condition:exon,
                                            flattenedfile) {
  .check_pkg("DEXSeq")
  .check_pkg("GenomicRanges")
  .check_pkg("IRanges")

  if (is.null(flattenedfile) || !file.exists(flattenedfile)) {
    stop("A valid flattened featureCounts GTF must be supplied.")
  }

  if (is.null(fc$counts) || is.null(fc$annotation)) {
    stop("featureCounts result must contain both 'counts' and 'annotation'.")
  }

  counts_df <- as.data.frame(fc$counts, check.names = FALSE)
  colnames(counts_df) <- rownames(sampleData)
  ann <- as.data.frame(fc$annotation, stringsAsFactors = FALSE)

  if (nrow(counts_df) != nrow(ann)) {
    stop("featureCounts annotation and count matrix have different numbers of rows.")
  }

  ord <- order(ann$GeneID, ann$Start, ann$End)
  counts_df <- counts_df[ord, , drop = FALSE]
  ann <- ann[ord, , drop = FALSE]

  id <- as.character(ann$GeneID)
  id <- substr(id, 1, 254)
  n <- ave(seq_along(id), id, FUN = seq_along)

  row_ids <- sprintf("%s:E%03d", id, as.numeric(n))
  if (anyDuplicated(row_ids)) {
    stop("Duplicate exon-bin IDs were produced after featureCounts truncation.")
  }
  rownames(counts_df) <- row_ids

  counts_mat <- as.matrix(counts_df)

  dcounts <- counts_mat[substr(rownames(counts_mat), 1, 1) != "_", , drop = FALSE]

  splitted <- strsplit(rownames(dcounts), ":")
  exons <- sapply(splitted, `[[`, 2)
  genesrle <- sapply(splitted, `[[`, 1)

  aggregates <- read.delim(flattenedfile, stringsAsFactors = FALSE, header = FALSE)
  colnames(aggregates) <- c(
    "chr", "source", "class", "start",
    "end", "score", "strand", "phase", "attr"
  )
  aggregates$strand <- gsub("\\.", "*", aggregates$strand)
  aggregates <- aggregates[aggregates$class == "exon", , drop = FALSE]

  aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
  aggregates$gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1", aggregates$attr)
  aggregates$gene_id <- substr(aggregates$gene_id, 1, 254)

  transcripts <- gsub(".*transcripts\\s(\\S+).*", "\\1", aggregates$attr)
  transcripts <- strsplit(transcripts, "\\+")
  exonids <- gsub(".*exon_number\\s(\\S+).*", "\\1", aggregates$attr)

  exoninfo <- GenomicRanges::GRanges(
    seqnames = as.character(aggregates$chr),
    ranges = IRanges::IRanges(start = aggregates$start, end = aggregates$end),
    strand = aggregates$strand
  )
  names(exoninfo) <- paste(aggregates$gene_id, exonids, sep = ":E")
  names(transcripts) <- names(exoninfo)

  if (!all(rownames(dcounts) %in% names(exoninfo))) {
    missing_ids <- setdiff(rownames(dcounts), names(exoninfo))
    stop(
      "Count rows do not correspond to flattened annotation file.\n",
      "Examples:\n",
      paste(utils::head(missing_ids, 10), collapse = "\n")
    )
  }

  matching <- match(rownames(dcounts), names(exoninfo))

  DEXSeq::DEXSeqDataSet(
    dcounts,
    sampleData,
    design,
    exons,
    genesrle,
    exoninfo[matching],
    transcripts[matching]
  )
}

.compare_methods <- function(overwrite = FALSE,
                             load_existing = TRUE) {
  outdir <- file.path(.get_results_dir(), "analysis/dexseq")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  policy_grid <- expand.grid(
    counting_method = c("summarizeOverlaps", "featureCounts"),
    count_multi_exon_overlaps = c(TRUE, FALSE),
    include_multimappers = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  run_table <- data.frame(
    counting_method = c(policy_grid$counting_method, "HTSeq"),
    count_multi_exon_overlaps = c(policy_grid$count_multi_exon_overlaps, TRUE),
    include_multimappers = c(policy_grid$include_multimappers, FALSE),
    stringsAsFactors = FALSE
  )

  mode_tags <- vapply(
    seq_len(nrow(run_table)),
    function(i) {
      .get_dexseq_mode_tag(
        counting_method = run_table$counting_method[i],
        count_multi_exon_overlaps = run_table$count_multi_exon_overlaps[i],
        include_multimappers = run_table$include_multimappers[i]
      )
    },
    character(1)
  )

  run_table$mode_tag <- mode_tags
  run_table$out_rds <- file.path(
    outdir,
    paste0("dexseq_results.H3.", run_table$mode_tag, ".rds")
  )
  run_table$out_xlsx <- file.path(
    outdir,
    paste0("dexseq_results.H3.", run_table$mode_tag, ".xlsx")
  )

  out <- vector("list", nrow(run_table))
  names(out) <- run_table$mode_tag

  for (i in seq_len(nrow(run_table))) {
    mode_tag <- run_table$mode_tag[i]
    out_rds <- run_table$out_rds[i]
    out_xlsx <- run_table$out_xlsx[i]

    if (file.exists(out_rds) && !overwrite) {
      if (isTRUE(load_existing)) {
        message("[DEXSeq H3 grid] loading existing: ", mode_tag)
        out[[mode_tag]] <- readRDS(out_rds)
      } else {
        message("[DEXSeq H3 grid] skipping existing: ", mode_tag)
        out[[mode_tag]] <- NULL
      }
      next
    }

    message("[DEXSeq H3 grid] running: ", mode_tag)

    out[[mode_tag]] <- run_dexseq(
      timepoints = "H3",
      counting_method = run_table$counting_method[i],
      count_multi_exon_overlaps = run_table$count_multi_exon_overlaps[i],
      include_multimappers = run_table$include_multimappers[i],
      force = overwrite,
      out_rds = out_rds,
      out_xlsx = out_xlsx
    )
  }

  out
}

.filter_dexseq_se_by_gene_universe <- function(se, gene_universe) {
  if (is.null(gene_universe)) {
    return(se)
  }

  gene_universe <- unique(.strip_ens_version(gene_universe))

  rd <- SummarizedExperiment::rowData(se)

  if ("gene_id" %in% names(rd)) {
    se_gene <- .strip_ens_version(rd$gene_id)
  } else if ("groupID" %in% names(rd)) {
    se_gene <- .strip_ens_version(rd$groupID)
  } else {
    se_gene <- .strip_ens_version(sub(":.*$", "", rownames(se)))
  }

  keep <- se_gene %in% gene_universe

  if (!any(keep)) {
    stop("DEXSeq filtered SE has zero rows after applying gene universe.")
  }

  se[keep, ]
}

.dexseq_count_dispatch <- function(method = c("summarizeOverlaps", "featureCounts", "HTSeq"),
                                   bam_files,
                                   flattened_annotation = NULL,
                                   sample_table,
                                   fc_gtf = NULL,
                                   htseq_count_files = NULL,
                                   htseq_gff = NULL,
                                   paired_end = TRUE,
                                   ignore_strand = TRUE,
                                   fragments = TRUE,
                                   include_multimappers = FALSE,
                                   count_multi_exon_overlaps = TRUE,
                                   yield_size = 1000000,
                                   threads = 1,
                                   BPPARAM = BiocParallel::SerialParam()) {
  method <- match.arg(method)

  if (method == "summarizeOverlaps") {
    return(.count_dexseq_summarizeOverlaps(
      bam_files = bam_files,
      flattened_annotation = flattened_annotation,
      sample_table = sample_table,
      paired_end = paired_end,
      ignore_strand = ignore_strand,
      fragments = fragments,
      include_multimappers = include_multimappers,
      count_multi_exon_overlaps = count_multi_exon_overlaps,
      yield_size = yield_size,
      BPPARAM = BPPARAM
    ))
  }

  if (method == "featureCounts") {
    if (is.null(fc_gtf)) {
      stop("fc_gtf must be supplied for featureCounts backend.")
    }

    return(.count_dexseq_featureCounts(
      bam_files = bam_files,
      fc_gtf = fc_gtf,
      paired_end = paired_end,
      ignore_strand = ignore_strand,
      count_multi_exon_overlaps = count_multi_exon_overlaps,
      include_multimappers = include_multimappers,
      threads = threads
    ))
  }

  if (method == "HTSeq") {
    if (is.null(htseq_count_files)) {
      stop("htseq_count_files must be supplied for HTSeq backend.")
    }
    if (is.null(htseq_gff)) {
      stop("htseq_gff must be supplied for HTSeq backend.")
    }

    return(.count_dexseq_HTSeq(
      count_files = htseq_count_files,
      sample_table = sample_table,
      flattenedfile = htseq_gff
    ))
  }

  stop("Unsupported counting method: ", method)
}

.run_one_dexseq_timepoint <- function(tp,
                                      coldata,
                                      flattened_annotation,
                                      counting_method = "summarizeOverlaps",
                                      contrast_var = "condition",
                                      ref_level = "C1",
                                      padj_cutoff = 0.10,
                                      paired_end = TRUE,
                                      ignore_strand = TRUE,
                                      fragments = TRUE,
                                      include_multimappers = FALSE,
                                      count_multi_exon_overlaps = TRUE,
                                      gene_universe_by_tp = NULL,
                                      yield_size = 1000000,
                                      BPPARAM,
                                      fc_gtf = NULL,
                                      htseq_gff = NULL,
                                      threads = 1,
                                      annot_df = NULL) {
  message("\n[DEXSeq] Starting analysis for: ", tp, " vs ", ref_level)

  cmp <- paste0(tp, "_vs_", ref_level)

  keep <- coldata[[contrast_var]] %in% c(ref_level, tp)
  csub <- coldata[keep, , drop = FALSE]

  sample_table <- data.frame(
    row.names = csub$names,
    sample = csub$names,
    condition = csub[[contrast_var]],
    libType = "paired-end",
    stringsAsFactors = FALSE
  )
  sample_table <- .ensure_dexseq_factors(sample_table, ref_level)

  bam_files <- NULL
  if ("bam" %in% colnames(csub)) {
    bam_files <- csub$bam
    names(bam_files) <- csub$names
  }

  htseq_count_files <- NULL
  if ("count_file" %in% colnames(csub)) {
    htseq_count_files <- csub$count_file
    names(htseq_count_files) <- csub$names
  }

  if (counting_method %in% c("summarizeOverlaps", "featureCounts") && is.null(bam_files)) {
    stop(
      "BAM-based counting selected ('", counting_method,
      "'), but no BAM files were found in coldata."
    )
  }

  if (counting_method == "HTSeq" && is.null(htseq_count_files)) {
    stop("HTSeq counting selected, but no HTSeq count files were found in coldata.")
  }

  if (counting_method == "HTSeq") {
    message("[DEXSeq] ", tp, ": using ", length(htseq_count_files), " HTSeq count files")
  } else {
    message("[DEXSeq] ", tp, ": using ", length(bam_files), " BAM files")
  }

  se <- NULL
  fc_res <- NULL

  timing <- system.time({
    message("[DEXSeq] ", tp, ": counting reads (", counting_method, ")")

    if (counting_method == "summarizeOverlaps") {
      se <- .dexseq_count_dispatch(
        method = counting_method,
        bam_files = bam_files,
        flattened_annotation = flattened_annotation,
        sample_table = sample_table,
        paired_end = paired_end,
        ignore_strand = ignore_strand,
        fragments = fragments,
        include_multimappers = include_multimappers,
        count_multi_exon_overlaps = count_multi_exon_overlaps,
        yield_size = yield_size,
        BPPARAM = BPPARAM,
        fc_gtf = fc_gtf,
        threads = threads
      )

      gene_universe <- NULL
      if (!is.null(gene_universe_by_tp)) {
        gene_universe <- gene_universe_by_tp[[tp]]
      }

      se <- .filter_dexseq_se_by_gene_universe(
        se = se,
        gene_universe = gene_universe
      )

      cd <- as.data.frame(SummarizedExperiment::colData(se))
      cd <- .ensure_dexseq_factors(cd, ref_level)
      SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(cd)

      message("[DEXSeq] ", tp, ": constructing DEXSeqDataSet from summarizeOverlaps")
      dxd <- suppressMessages(
        suppressWarnings(
          DEXSeq::DEXSeqDataSetFromSE(
            se,
            design = ~ sample + exon + condition:exon
          )
        )
      )
    } else if (counting_method == "featureCounts") {
      fc_res <- .dexseq_count_dispatch(
        method = counting_method,
        bam_files = bam_files,
        flattened_annotation = flattened_annotation,
        sample_table = sample_table,
        paired_end = paired_end,
        ignore_strand = ignore_strand,
        fragments = fragments,
        include_multimappers = include_multimappers,
        count_multi_exon_overlaps = count_multi_exon_overlaps,
        yield_size = yield_size,
        BPPARAM = BPPARAM,
        fc_gtf = fc_gtf,
        threads = threads
      )

      message("[DEXSeq] ", tp, ": constructing DEXSeqDataSet from featureCounts")
      dxd <- suppressWarnings(
        .DEXSeqDataSetFromFeatureCounts(
          fc = fc_res,
          sampleData = sample_table,
          design = ~ sample + exon + condition:exon,
          flattenedfile = fc_gtf
        )
      )

      se <- NULL
    } else if (counting_method == "HTSeq") {
      dxd <- .dexseq_count_dispatch(
        method = counting_method,
        bam_files = bam_files,
        flattened_annotation = flattened_annotation,
        sample_table = sample_table,
        fc_gtf = fc_gtf,
        htseq_count_files = htseq_count_files,
        htseq_gff = htseq_gff,
        paired_end = paired_end,
        ignore_strand = ignore_strand,
        fragments = fragments,
        include_multimappers = include_multimappers,
        count_multi_exon_overlaps = count_multi_exon_overlaps,
        yield_size = yield_size,
        threads = threads,
        BPPARAM = BPPARAM
      )

      message("[DEXSeq] ", tp, ": using DEXSeqDataSetFromHTSeq input")
      se <- NULL
      fc_res <- NULL
    } else {
      stop("Unsupported counting method: ", counting_method)
    }

    message(
      "[DEXSeq] ", tp, ": modeling ",
      nrow(dxd), " exon bins across ",
      ncol(dxd), " DEXSeq model columns from ",
      nrow(sample_table), " biological samples"
    )

    message("[DEXSeq] ", tp, ": estimating size factors")
    dxd <- DEXSeq::estimateSizeFactors(dxd)

    message("[DEXSeq] ", tp, ": estimating dispersions")
    dxd <- DEXSeq::estimateDispersions(dxd)

    message("[DEXSeq] ", tp, ": testing for differential exon usage")
    dxd <- DEXSeq::testForDEU(dxd)

    message("[DEXSeq] ", tp, ": estimating exon fold changes")
    dxd <- DEXSeq::estimateExonFoldChanges(
      dxd,
      fitExpToVar = "condition"
    )

    message("[DEXSeq] ", tp, ": collecting results")
    dxr <- DEXSeq::DEXSeqResults(dxd)
  })

  message("[DEXSeq] ", tp, ": runtime = ", round(timing["elapsed"], 2), " sec")

  exon_full <- .standardize_dexseq_exon_res(
    res_df = dxr,
    timepoint = tp,
    contrast = cmp,
    annot_df = annot_df
  )

  gene_q <- DEXSeq::perGeneQValue(dxr)
  gene_tab <- data.frame(
    groupID = names(gene_q),
    padj = as.numeric(gene_q),
    stringsAsFactors = FALSE
  )

  gene_full <- .standardize_dexseq_gene_res(
    res_df = gene_tab,
    timepoint = tp,
    contrast = cmp,
    annot_df = annot_df
  )

  exon_full$significant <- !is.na(exon_full$padj) & exon_full$padj <= padj_cutoff
  gene_full$significant <- !is.na(gene_full$padj) & gene_full$padj <= padj_cutoff

  exon_sig <- exon_full[exon_full$significant, , drop = FALSE]
  gene_sig <- gene_full[gene_full$significant, , drop = FALSE]

  if ("padj" %in% colnames(gene_sig)) {
    gene_sig <- gene_sig[order(gene_sig$padj), , drop = FALSE]
  }
  if ("padj" %in% colnames(exon_sig)) {
    exon_sig <- exon_sig[order(exon_sig$padj), , drop = FALSE]
  }

  list(
    dxd = dxd,
    dxr = dxr,
    se = se,
    counts = if (!is.null(se)) SummarizedExperiment::assay(se) else DEXSeq::counts(dxd),
    fc = fc_res,
    results = list(
      gene_full = gene_full,
      gene_sig = gene_sig,
      exon_full = exon_full,
      exon_sig = exon_sig
    ),
    filtered_genes = unique(.strip_ens_version(gene_full$gene_id_full))
  )
}

# -------------------------
# Public API
# -------------------------

run_dexseq <- function(samples_tsv = NULL,
                       env_file = NULL,
                       resources_dir = NULL,
                       timepoints = c("H1", "H3", "H24"),
                       contrast_var = "condition",
                       ref_level = "C1",
                       padj_cutoff = 0.10,
                       counting_method = "summarizeOverlaps",
                       paired_end = TRUE,
                       ignore_strand = TRUE,
                       fragments = TRUE,
                       include_multimappers = TRUE,
                       count_multi_exon_overlaps = TRUE,
                       gene_universe_by_tp = NULL,
                       annot_df = NULL,
                       workers = NULL,
                       yield_size = 1000000,
                       force = FALSE,
                       out_rds = NULL,
                       out_xlsx = NULL) {
  .check_pkg("DEXSeq")
  .check_pkg("DESeq2")
  .check_pkg("GenomicFeatures")
  .check_pkg("GenomicAlignments")
  .check_pkg("Rsamtools")
  .check_pkg("BiocParallel")
  .check_pkg("openxlsx")

  mode_tag <- .get_dexseq_mode_tag(
    counting_method = counting_method,
    count_multi_exon_overlaps = count_multi_exon_overlaps,
    include_multimappers = include_multimappers
  )

  outdir <- file.path(.get_results_dir(), "analysis/dexseq")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(out_rds)) {
    out_rds <- file.path(outdir, paste0("dexseq_results.", ref_level, ".", mode_tag, ".rds"))
  }
  if (is.null(out_xlsx)) {
    out_xlsx <- file.path(outdir, paste0("dexseq_results.", ref_level, ".", mode_tag, ".xlsx"))
  }

  if (file.exists(out_rds) && !force) {
    message("Loading cached DEXSeq results: ", out_rds)
    return(readRDS(out_rds))
  }

  paths <- .get_runtime_paths()

  samples_tsv <- samples_tsv %||% paths$samples_tsv
  env_file <- env_file %||% paths$env_file
  resources_dir <- resources_dir %||% paths$resources_dir

  if (is.null(workers)) {
    threads_env <- .resolve_threads_from_env(env_file)
    workers <- max(1, min(4, floor(threads_env / 2)))
  }
  if (is.null(workers) || is.na(workers)) workers <- 1

  if (is.null(annot_df)) {
    annot_df <- .get_annot()
  }

  if (counting_method == "HTSeq") {
    coldata <- .read_dexseq_htseq_coldata(
      samples_tsv = samples_tsv,
      results_dir = .get_results_dir(),
      contrast_var = contrast_var,
      sample_col = "sample_id"
    )
  } else {
    coldata <- .read_dexseq_bam_coldata(
      samples_tsv = samples_tsv,
      results_dir = .get_results_dir(),
      contrast_var = contrast_var,
      sample_col = "sample_id"
    )
  }

  flattened_annotation <- NULL
  fc_gtf <- NULL
  htseq_gff <- NULL

  if (counting_method == "summarizeOverlaps") {
    flattened_annotation <- .get_dexseq_flattened_annotation(
      env_file = env_file,
      resources_dir = resources_dir,
      overwrite = FALSE
    )
  }

  if (counting_method == "featureCounts") {
    fc_gtf <- .get_dexseq_featurecounts_gtf(
      env_file = env_file,
      resources_dir = resources_dir,
      overwrite = FALSE
    )
  }

  if (counting_method == "HTSeq") {
    htseq_gff <- .get_dexseq_htseq_gff(
      env_file = env_file,
      resources_dir = resources_dir,
      overwrite = FALSE
    )
  }

  BPPARAM <- if (.Platform$OS.type == "windows") {
    BiocParallel::SnowParam(workers = workers, progressbar = TRUE)
  } else {
    BiocParallel::MulticoreParam(workers = workers, progressbar = TRUE)
  }

  by_tp <- setNames(lapply(timepoints, function(tp) {
    .run_one_dexseq_timepoint(
      tp = tp,
      coldata = coldata,
      flattened_annotation = flattened_annotation,
      counting_method = counting_method,
      contrast_var = contrast_var,
      ref_level = ref_level,
      padj_cutoff = padj_cutoff,
      paired_end = paired_end,
      ignore_strand = ignore_strand,
      fragments = fragments,
      include_multimappers = include_multimappers,
      count_multi_exon_overlaps = count_multi_exon_overlaps,
      gene_universe_by_tp = gene_universe_by_tp,
      yield_size = yield_size,
      BPPARAM = BPPARAM,
      fc_gtf = fc_gtf,
      htseq_gff = htseq_gff,
      threads = workers,
      annot_df = annot_df
    )
  }), timepoints)

  gene_full_list <- lapply(by_tp, function(x) x$results$gene_full)
  gene_sig_list <- lapply(by_tp, function(x) x$results$gene_sig)
  exon_full_list <- lapply(by_tp, function(x) x$results$exon_full)
  exon_sig_list <- lapply(by_tp, function(x) x$results$exon_sig)

  summary_tbl <- .make_dexseq_summary(
    gene_sig_list = gene_sig_list,
    exon_sig_list = exon_sig_list,
    padj_cutoff = padj_cutoff
  )

  .write_dexseq_xlsx(
    gene_sig_list = gene_sig_list,
    exon_sig_list = exon_sig_list,
    summary_tbl = summary_tbl,
    out_xlsx = out_xlsx
  )

  out <- list(
    meta = list(
      project_root = PROJECT_ROOT,
      samples_tsv = samples_tsv,
      env_file = env_file,
      resources_dir = resources_dir,
      contrast_var = contrast_var,
      ref_level = ref_level,
      timepoints = timepoints,
      padj_cutoff = padj_cutoff,
      counting_method = counting_method,
      paired_end = paired_end,
      ignore_strand = ignore_strand,
      fragments = fragments,
      include_multimappers = include_multimappers,
      count_multi_exon_overlaps = count_multi_exon_overlaps,
      workers = workers,
      yield_size = yield_size,
      gene_universe_by_tp = gene_universe_by_tp,
      gene_universe_n_by_tp = if (!is.null(gene_universe_by_tp)) {
        vapply(gene_universe_by_tp, function(x) length(unique(.strip_ens_version(x))), numeric(1))
      } else {
        NULL
      }
    ),
    by_tp = by_tp,
    results = list(
      gene_full_by_tp = gene_full_list,
      gene_sig_by_tp = gene_sig_list,
      exon_full_by_tp = exon_full_list,
      exon_sig_by_tp = exon_sig_list,
      gene_full_all = .rbind_fill(gene_full_list),
      gene_sig_all = .rbind_fill(gene_sig_list),
      exon_full_all = .rbind_fill(exon_full_list),
      exon_sig_all = .rbind_fill(exon_sig_list),
      sig_all = list(
        genes = .rbind_fill(gene_sig_list),
        exons = .rbind_fill(exon_sig_list)
      )
    ),
    summary = summary_tbl,
    paths = list(
      rds = out_rds,
      xlsx = out_xlsx
    )
  )

  dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
  saveRDS(out, out_rds)

  out
}
