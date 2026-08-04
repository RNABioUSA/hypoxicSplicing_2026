#!/usr/bin/env Rscript

parse_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    if (!startsWith(key, "--") || i == length(x)) {
      stop("Arguments must be supplied as --name value pairs; problem near: ", key)
    }
    out[[sub("^--", "", key)]] <- x[[i + 1L]]
    i <- i + 2L
  }
  out
}

`%||%` <- function(x, y) if (is.null(x)) y else x

args <- parse_args(commandArgs(trailingOnly = TRUE))
required_args <- c("metadata", "quant-root", "gtf", "out")
missing_args <- required_args[!vapply(required_args, function(z) nzchar(args[[z]] %||% ""), logical(1))]
if (length(missing_args)) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

metadata_path <- normalizePath(args[["metadata"]], mustWork = TRUE)
quant_root <- normalizePath(args[["quant-root"]], mustWork = TRUE)
gtf_path <- normalizePath(args[["gtf"]], mustWork = TRUE)
out_root <- normalizePath(args[["out"]], mustWork = FALSE)
n_perm <- as.integer(args[["n-perm"]] %||% 20L)
seed <- as.integer(args[["seed"]] %||% 20260727L)
alpha <- as.numeric(args[["alpha"]] %||% 0.10)
workers <- as.integer(args[["workers"]] %||% 1L)

if (!is.finite(n_perm) || n_perm < 1L) stop("--n-perm must be a positive integer.")
if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) stop("--alpha must lie between 0 and 1.")
if (!is.finite(workers) || workers < 1L) workers <- 1L

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "cache"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "permutations"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "suppa"), recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "tximport", "DESeq2", "DRIMSeq", "stageR", "DEXSeq",
  "rtracklayer", "GenomicFeatures", "GenomicAlignments",
  "Rsamtools", "SummarizedExperiment", "S4Vectors", "BiocParallel"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages)) {
  stop(
    "Missing R/Bioconductor packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x, file = path, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
  invisible(path)
}

write_tsv_gz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(
    x, file = con, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
  invisible(path)
}

strip_version <- function(x) sub("\\.[0-9]+$", "", as.character(x))

subset_txi <- function(txi, sample_ids) {
  out <- txi
  for (nm in intersect(c("abundance", "counts", "length"), names(out))) {
    out[[nm]] <- out[[nm]][, sample_ids, drop = FALSE]
  }
  out
}

write_suppa_matrix <- function(mat, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(paste(colnames(mat), collapse = "\t"), con = con)
  utils::write.table(
    data.frame(transcript_id = rownames(mat), mat, check.names = FALSE),
    file = con, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE
  )
  invisible(path)
}

make_population_allocations <- function(ids) {
  if (length(ids) != 5L) {
    stop("Expected five 0h normoxia donors per population; found ", length(ids), ".")
  }
  a_sets <- combn(sort(ids), 2L, simplify = FALSE)
  rows <- list()
  k <- 0L
  for (a in a_sets) {
    remaining <- setdiff(ids, a)
    b_sets <- combn(sort(remaining), 2L, simplify = FALSE)
    for (b in b_sets) {
      k <- k + 1L
      rows[[k]] <- list(A = sort(a), B = sort(b))
    }
  }
  rows
}

make_permutations <- function(metadata, n, random_seed) {
  by_population <- split(metadata$sample_id, metadata$population)
  if (length(by_population) != 2L) {
    stop("The balanced design expects exactly two populations.")
  }
  pop_names <- sort(names(by_population))
  alloc_1 <- make_population_allocations(by_population[[pop_names[[1]]]])
  alloc_2 <- make_population_allocations(by_population[[pop_names[[2]]]])

  candidates <- vector("list", length(alloc_1) * length(alloc_2))
  canonical_keys <- character(length(candidates))
  k <- 0L
  for (i in seq_along(alloc_1)) {
    for (j in seq_along(alloc_2)) {
      k <- k + 1L
      group_a <- sort(c(alloc_1[[i]]$A, alloc_2[[j]]$A))
      group_b <- sort(c(alloc_1[[i]]$B, alloc_2[[j]]$B))
      key_ab <- paste0("A=", paste(group_a, collapse = ","), "|B=", paste(group_b, collapse = ","))
      key_ba <- paste0("A=", paste(group_b, collapse = ","), "|B=", paste(group_a, collapse = ","))
      canonical_keys[[k]] <- min(key_ab, key_ba)
      candidates[[k]] <- list(A = group_a, B = group_b)
    }
  }

  candidates <- candidates[!duplicated(canonical_keys)]
  if (n > length(candidates)) {
    stop("Requested ", n, " permutations, but only ", length(candidates), " unique balanced splits exist.")
  }

  set.seed(random_seed)
  selected <- sample(seq_along(candidates), n, replace = FALSE)
  candidates <- candidates[selected]

  membership <- do.call(rbind, lapply(seq_along(candidates), function(i) {
    z <- candidates[[i]]
    data.frame(
      perm_id = sprintf("perm_%03d", i),
      sample_id = c(z$A, z$B),
      group = rep(c("A", "B"), each = 4L),
      stringsAsFactors = FALSE
    )
  }))
  membership <- merge(
    membership,
    metadata[, c("sample_id", "population", "individual", "gsm", "run_accession")],
    by = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )
  membership <- membership[
    order(membership$perm_id, membership$group, membership$population, membership$sample_id),
    c("perm_id", "sample_id", "group", "population", "individual", "gsm", "run_accession")
  ]
  rownames(membership) <- NULL
  membership
}

metadata <- utils::read.delim(
  metadata_path, stringsAsFactors = FALSE, check.names = FALSE
)
required_columns <- c("sample_id", "population", "individual", "gsm", "run_accession", "stress", "time")
missing_columns <- setdiff(required_columns, colnames(metadata))
if (length(missing_columns)) {
  stop("Metadata is missing columns: ", paste(missing_columns, collapse = ", "))
}
if (anyDuplicated(metadata$sample_id)) stop("sample_id values must be unique.")
if (nrow(metadata) != 10L) stop("Expected ten 0h normoxia samples; found ", nrow(metadata), ".")
if (any(metadata$stress != "normoxia") || any(metadata$time != "0h")) {
  stop("This workflow accepts only untreated 0h normoxia samples.")
}

salmon_files <- file.path(quant_root, "salmon", metadata$sample_id, "quant.sf")
names(salmon_files) <- metadata$sample_id
bam_files <- file.path(
  quant_root, "star", metadata$sample_id, "Aligned.sortedByCoord.out.bam"
)
names(bam_files) <- metadata$sample_id

missing_inputs <- c(salmon_files[!file.exists(salmon_files)], bam_files[!file.exists(bam_files)])
if (length(missing_inputs)) {
  stop("Quantification/alignment files are missing:\n", paste(missing_inputs, collapse = "\n"))
}

membership <- make_permutations(metadata, n_perm, seed)
write_tsv(membership, file.path(out_root, "permutation_membership.tsv"))

message("[setup] importing GTF and creating transcript-to-gene map")
gtf <- rtracklayer::import(gtf_path)
tx_ids <- as.character(S4Vectors::mcols(gtf)$transcript_id)
gene_ids <- as.character(S4Vectors::mcols(gtf)$gene_id)
keep_tx <- !is.na(tx_ids) & nzchar(tx_ids) & !is.na(gene_ids) & nzchar(gene_ids)
tx2gene <- unique(data.frame(
  TXNAME = tx_ids[keep_tx],
  GENEID = gene_ids[keep_tx],
  stringsAsFactors = FALSE
))
rm(gtf)

message("[setup] importing Salmon estimates")
txi_gene <- tximport::tximport(
  salmon_files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = FALSE
)
txi_tx <- tximport::tximport(
  salmon_files,
  type = "salmon",
  txOut = TRUE,
  countsFromAbundance = "dtuScaledTPM"
)
saveRDS(txi_gene, file.path(out_root, "cache", "tximport_gene.rds"))
saveRDS(txi_tx, file.path(out_root, "cache", "tximport_transcript.rds"))
write_suppa_matrix(
  txi_tx$abundance,
  file.path(out_root, "suppa", "all_samples_iso_tpm.txt")
)

tx_match <- match(rownames(txi_tx$counts), tx2gene$TXNAME)
if (anyNA(tx_match)) {
  stop(sum(is.na(tx_match)), " quantified transcripts could not be mapped to a gene in the supplied GTF.")
}
drim_counts_all <- data.frame(
  gene_id = tx2gene$GENEID[tx_match],
  feature_id = rownames(txi_tx$counts),
  txi_tx$counts,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

message("[setup] building/loading DEXSeq exon-bin counts")
dexseq_cache <- file.path(out_root, "cache", "dexseq_exon_bin_counts.rds")
if (file.exists(dexseq_cache)) {
  dexseq_se_all <- readRDS(dexseq_cache)
} else {
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_path, format = "gtf")
  flattened <- GenomicFeatures::exonicParts(
    txdb, linked.to.single.gene.only = TRUE
  )
  flat_gene <- as.character(S4Vectors::mcols(flattened)$gene_id)
  flat_part <- as.integer(S4Vectors::mcols(flattened)$exonic_part)
  names(flattened) <- sprintf("%s:E%03d", flat_gene, flat_part)

  bam_list <- Rsamtools::BamFileList(lapply(
    bam_files,
    function(z) Rsamtools::BamFile(z, yieldSize = 1000000L)
  ))
  bp <- if (.Platform$OS.type == "windows") {
    BiocParallel::SnowParam(workers = workers, progressbar = TRUE)
  } else {
    BiocParallel::MulticoreParam(workers = workers, progressbar = TRUE)
  }
  dexseq_se_all <- GenomicAlignments::summarizeOverlaps(
    features = flattened,
    reads = bam_list,
    singleEnd = FALSE,
    fragments = TRUE,
    ignore.strand = TRUE,
    inter.feature = FALSE,
    param = Rsamtools::ScanBamParam(tag = "NH"),
    BPPARAM = bp
  )
  colnames(dexseq_se_all) <- names(bam_files)
  saveRDS(dexseq_se_all, dexseq_cache)
}

bp <- if (.Platform$OS.type == "windows") {
  BiocParallel::SnowParam(workers = workers, progressbar = FALSE)
} else {
  BiocParallel::MulticoreParam(workers = workers, progressbar = FALSE)
}

summary_rows <- list()
significant_gene_rows <- list()
summary_i <- 0L
sig_i <- 0L

add_summary <- function(perm_id, tool, level, n_tested, n_significant, threshold, status = "ok") {
  summary_i <<- summary_i + 1L
  summary_rows[[summary_i]] <<- data.frame(
    perm_id = perm_id,
    tool = tool,
    level = level,
    n_tested = as.integer(n_tested),
    n_significant = as.integer(n_significant),
    fraction_significant = if (is.finite(n_tested) && n_tested > 0) n_significant / n_tested else NA_real_,
    threshold = threshold,
    status = status,
    stringsAsFactors = FALSE
  )
}

add_sig_genes <- function(perm_id, tool, ids) {
  ids <- sort(unique(strip_version(ids)))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (!length(ids)) return(invisible(NULL))
  sig_i <<- sig_i + 1L
  significant_gene_rows[[sig_i]] <<- data.frame(
    perm_id = perm_id,
    tool = tool,
    gene_id = ids,
    stringsAsFactors = FALSE
  )
  invisible(NULL)
}

perm_ids <- unique(membership$perm_id)
for (perm_index in seq_along(perm_ids)) {
  perm_id <- perm_ids[[perm_index]]
  message("\n[permutation] ", perm_id, " (", perm_index, "/", length(perm_ids), ")")
  pmeta <- membership[membership$perm_id == perm_id, , drop = FALSE]
  sample_ids <- pmeta$sample_id
  condition <- factor(pmeta$group, levels = c("A", "B"))
  sample_data <- data.frame(
    row.names = sample_ids,
    condition = condition,
    stringsAsFactors = TRUE
  )
  perm_out <- file.path(out_root, "permutations", perm_id)
  dir.create(perm_out, recursive = TRUE, showWarnings = FALSE)

  tryCatch({
    message("[DESeq2] ", perm_id)
    dds <- DESeq2::DESeqDataSetFromTximport(
      subset_txi(txi_gene, sample_ids),
      colData = sample_data,
      design = ~condition
    )
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    deseq_res <- DESeq2::results(
      dds,
      contrast = c("condition", "B", "A"),
      alpha = alpha
    )
    deseq_df <- as.data.frame(deseq_res)
    deseq_df$gene_id <- rownames(deseq_df)
    rownames(deseq_df) <- NULL
    deseq_df <- deseq_df[, c("gene_id", setdiff(colnames(deseq_df), "gene_id"))]
    write_tsv_gz(deseq_df, file.path(perm_out, "deseq2_gene_results.tsv.gz"))

    tested <- !is.na(deseq_df$padj)
    significant <- tested & deseq_df$padj < alpha
    add_summary(
      perm_id, "DESeq2", "gene",
      sum(tested), sum(significant),
      sprintf("BH FDR < %.3g", alpha)
    )
    add_sig_genes(perm_id, "DESeq2", deseq_df$gene_id[significant])
  }, error = function(e) {
    warning("[DESeq2] ", perm_id, " failed: ", conditionMessage(e))
    add_summary(perm_id, "DESeq2", "gene", NA, NA, "error", conditionMessage(e))
  })

  tryCatch({
    message("[DRIMSeq] ", perm_id)
    count_cols <- match(sample_ids, colnames(drim_counts_all))
    counts_df <- drim_counts_all[, c(1L, 2L, count_cols), drop = FALSE]
    keep_nonzero <- rowSums(counts_df[, sample_ids, drop = FALSE]) > 0
    counts_df <- counts_df[keep_nonzero, , drop = FALSE]
    samples_df <- data.frame(
      sample_id = sample_ids,
      condition = condition,
      stringsAsFactors = FALSE
    )
    samples_df$condition <- factor(samples_df$condition, levels = c("A", "B"))

    d <- DRIMSeq::dmDSdata(counts = counts_df, samples = samples_df)
    d <- DRIMSeq::dmFilter(
      d,
      min_samps_feature_expr = 3,
      min_samps_gene_expr = 6,
      min_feature_expr = 10,
      min_gene_expr = 10
    )
    design <- stats::model.matrix(~condition, data = DRIMSeq::samples(d))
    set.seed(seed + perm_index)
    d <- DRIMSeq::dmPrecision(d, design = design, BPPARAM = bp)
    d <- DRIMSeq::dmFit(d, design = design, BPPARAM = bp)
    d <- DRIMSeq::dmTest(d, coef = "conditionB")

    gene_res <- as.data.frame(DRIMSeq::results(d))
    tx_res <- as.data.frame(DRIMSeq::results(d, level = "feature"))
    p_screen <- gene_res$pvalue
    p_screen[is.na(p_screen)] <- 1
    names(p_screen) <- gene_res$gene_id
    p_confirmation <- matrix(tx_res$pvalue, ncol = 1L)
    p_confirmation[is.na(p_confirmation)] <- 1
    rownames(p_confirmation) <- tx_res$feature_id

    stage <- stageR::stageRTx(
      pScreen = p_screen,
      pConfirmation = p_confirmation,
      pScreenAdjusted = FALSE,
      tx2gene = tx_res[, c("feature_id", "gene_id")]
    )
    stage <- stageR::stageWiseAdjustment(stage, method = "dtu", alpha = alpha)
    sig_gene_ids <- rownames(as.data.frame(stageR::getSignificantGenes(stage)))
    sig_tx_ids <- rownames(as.data.frame(stageR::getSignificantTx(stage)))

    gene_res$stageR_significant <- gene_res$gene_id %in% sig_gene_ids
    tx_res$stageR_significant <- tx_res$feature_id %in% sig_tx_ids
    write_tsv_gz(gene_res, file.path(perm_out, "drimseq_gene_results.tsv.gz"))
    write_tsv_gz(tx_res, file.path(perm_out, "drimseq_transcript_results.tsv.gz"))

    add_summary(
      perm_id, "DRIMSeq", "gene",
      sum(!is.na(gene_res$pvalue)), length(unique(sig_gene_ids)),
      sprintf("stageR OFDR %.3g", alpha)
    )
    add_summary(
      perm_id, "DRIMSeq", "transcript",
      sum(!is.na(tx_res$pvalue)), length(unique(sig_tx_ids)),
      sprintf("stageR OFDR %.3g", alpha)
    )
    add_sig_genes(perm_id, "DRIMSeq", sig_gene_ids)
  }, error = function(e) {
    warning("[DRIMSeq] ", perm_id, " failed: ", conditionMessage(e))
    add_summary(perm_id, "DRIMSeq", "gene", NA, NA, "error", conditionMessage(e))
    add_summary(perm_id, "DRIMSeq", "transcript", NA, NA, "error", conditionMessage(e))
  })

  tryCatch({
    message("[DEXSeq] ", perm_id)
    se <- dexseq_se_all[, sample_ids, drop = FALSE]
    dex_sample_data <- data.frame(
      row.names = sample_ids,
      sample = factor(sample_ids),
      condition = condition,
      libType = "paired-end",
      stringsAsFactors = TRUE
    )
    SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(dex_sample_data)
    dxd <- suppressMessages(suppressWarnings(
      DEXSeq::DEXSeqDataSetFromSE(
        se,
        design = ~sample + exon + condition:exon
      )
    ))
    dxd <- DEXSeq::estimateSizeFactors(dxd)
    dxd <- DEXSeq::estimateDispersions(dxd)
    dxd <- DEXSeq::testForDEU(dxd)
    dxr <- DEXSeq::DEXSeqResults(dxd)
    dex_df <- as.data.frame(dxr)
    write_tsv_gz(dex_df, file.path(perm_out, "dexseq_exon_results.tsv.gz"))

    gene_q <- DEXSeq::perGeneQValue(dxr)
    gene_df <- data.frame(
      gene_id = names(gene_q),
      padj = as.numeric(gene_q),
      stringsAsFactors = FALSE
    )
    write_tsv_gz(gene_df, file.path(perm_out, "dexseq_gene_results.tsv.gz"))

    tested_exons <- !is.na(dex_df$padj)
    sig_exons <- tested_exons & dex_df$padj <= alpha
    tested_genes <- !is.na(gene_df$padj)
    sig_genes <- tested_genes & gene_df$padj <= alpha
    add_summary(
      perm_id, "DEXSeq", "exon_bin",
      sum(tested_exons), sum(sig_exons),
      sprintf("BH FDR <= %.3g", alpha)
    )
    add_summary(
      perm_id, "DEXSeq", "gene",
      sum(tested_genes), sum(sig_genes),
      sprintf("perGeneQValue <= %.3g", alpha)
    )
    add_sig_genes(perm_id, "DEXSeq", gene_df$gene_id[sig_genes])
  }, error = function(e) {
    warning("[DEXSeq] ", perm_id, " failed: ", conditionMessage(e))
    add_summary(perm_id, "DEXSeq", "exon_bin", NA, NA, "error", conditionMessage(e))
    add_summary(perm_id, "DEXSeq", "gene", NA, NA, "error", conditionMessage(e))
  })
}

summary_df <- do.call(rbind, summary_rows)
write_tsv(summary_df, file.path(out_root, "statistics_summary.tsv"))

if (length(significant_gene_rows)) {
  sig_gene_df <- do.call(rbind, significant_gene_rows)
} else {
  sig_gene_df <- data.frame(
    perm_id = character(), tool = character(), gene_id = character(),
    stringsAsFactors = FALSE
  )
}
write_tsv(sig_gene_df, file.path(out_root, "significant_gene_ids.tsv"))

session_path <- file.path(out_root, "sessionInfo.txt")
session_con <- file(session_path, open = "wt")
writeLines(capture.output(sessionInfo()), con = session_con)
close(session_con)

message("\n[complete] R statistics written to: ", out_root)
