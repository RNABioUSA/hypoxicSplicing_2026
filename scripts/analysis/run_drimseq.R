PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
if (!nzchar(PROJECT_ROOT)) stop("PROJECT_ROOT env var not set.")

HELPERS_FILE <- file.path(PROJECT_ROOT, "scripts/utils/helpers.R")
if (!file.exists(HELPERS_FILE)) stop("helpers.R not found at: ", HELPERS_FILE)
source(HELPERS_FILE)

# -------------------------
# Internal helpers
# -------------------------

.standardize_drimseq_gene_res <- function(res_df,
                                          timepoint,
                                          contrast,
                                          tool = "drimseq",
                                          annot_df = NULL) {
  out <- res_df

  out$gene_id_full <- out$gene_id
  out$gene_id <- .strip_ens_version(out$gene_id)
  out$ensgene <- out$gene_id

  out$timepoint <- timepoint
  out$contrast <- contrast
  out$tool <- tool
  out$comparison <- contrast

  if ("adj_pvalue" %in% colnames(out)) {
    out$padj <- out$adj_pvalue
  } else {
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

.standardize_drimseq_tx_res <- function(res_df,
                                        timepoint,
                                        contrast,
                                        tool = "drimseq",
                                        annot_df = NULL) {
  out <- res_df

  out$gene_id_full <- out$gene_id
  out$gene_id <- .strip_ens_version(out$gene_id)
  out$ensgene <- out$gene_id

  out$enstx <- out$feature_id
  out$timepoint <- timepoint
  out$contrast <- contrast
  out$tool <- tool
  out$comparison <- contrast

  if ("adj_pvalue" %in% colnames(out)) {
    out$padj <- out$adj_pvalue
  } else {
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

.make_drimseq_summary <- function(gene_sig_list,
                                  tx_sig_list,
                                  padj_cutoff = 0.10) {
  tps <- names(gene_sig_list)
  do.call(rbind, lapply(tps, function(tp) {
    data.frame(
      timepoint = tp,
      padj_cutoff = padj_cutoff,
      n_sig_genes = nrow(gene_sig_list[[tp]]),
      n_sig_tx = nrow(tx_sig_list[[tp]]),
      stringsAsFactors = FALSE
    )
  }))
}

.write_drimseq_xlsx <- function(gene_sig_list,
                                tx_sig_list,
                                summary_tbl,
                                out_xlsx) {
  .check_pkg("openxlsx")

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Summary")
  openxlsx::writeData(wb, "Summary", summary_tbl)

  for (tp in names(gene_sig_list)) {
    gene_sheet <- paste0(tp, "_genes")
    tx_sheet <- paste0(tp, "_tx")

    openxlsx::addWorksheet(wb, gene_sheet)
    if (nrow(gene_sig_list[[tp]]) == 0) {
      openxlsx::writeData(wb, gene_sheet, data.frame(note = "No significant genes at this cutoff."))
    } else {
      openxlsx::writeData(wb, gene_sheet, gene_sig_list[[tp]])
    }

    openxlsx::addWorksheet(wb, tx_sheet)
    if (nrow(tx_sig_list[[tp]]) == 0) {
      openxlsx::writeData(wb, tx_sheet, data.frame(note = "No significant transcripts at this cutoff."))
    } else {
      openxlsx::writeData(wb, tx_sheet, tx_sig_list[[tp]])
    }
  }

  dir.create(dirname(out_xlsx), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
  invisible(TRUE)
}

.run_one_drimseq_timepoint <- function(tp,
                                       coldata,
                                       counts_all,
                                       contrast_var = "condition",
                                       ref_level = "C1",
                                       filter_params,
                                       padj_cutoff = 0.10,
                                       drop_single_tx_genes = FALSE,
                                       BPPARAM,
                                       annot_df = NULL) {
  message("\n[DRIMSeq] Starting analysis for: ", tp, " vs ", ref_level)

  cmp <- paste0(tp, "_vs_", ref_level)

  keep <- coldata[[contrast_var]] %in% c(ref_level, tp)
  csub <- coldata[keep, , drop = FALSE]

  samples_df <- data.frame(
    sample_id = csub$names,
    condition = csub[[contrast_var]],
    stringsAsFactors = FALSE
  )

  samples_df$condition <- factor(samples_df$condition)
  if (!(ref_level %in% levels(samples_df$condition))) {
    stop("Reference level not found in ", contrast_var, ": ", ref_level)
  }
  samples_df$condition <- stats::relevel(samples_df$condition, ref = ref_level)

  sample_ids <- csub$names

  count_idx <- match(sample_ids, colnames(counts_all))
  if (anyNA(count_idx)) {
    stop("Not all samples for ", tp, " were found in imported count matrix.")
  }
  count_cols <- colnames(counts_all)[count_idx]

  counts_df <- counts_all[, c("gene_id", "feature_id", count_cols), drop = FALSE]

  # Drop features with zero counts in this subset
  keep_rows <- rowSums(counts_df[, count_cols, drop = FALSE]) > 0
  counts_df <- counts_df[keep_rows, , drop = FALSE]

  if (drop_single_tx_genes) {
    tx_per_gene <- table(counts_df$gene_id)
    multi_tx_genes <- names(tx_per_gene[tx_per_gene > 1])
    counts_df <- counts_df[counts_df$gene_id %in% multi_tx_genes, , drop = FALSE]
  }

  d <- DRIMSeq::dmDSdata(counts = counts_df, samples = samples_df)
  d <- do.call(DRIMSeq::dmFilter, c(list(d), filter_params))

  design <- stats::model.matrix(~condition, data = DRIMSeq::samples(d))
  coef_name <- paste0("condition", tp)

  if (!(coef_name %in% colnames(design))) {
    stop(
      "Coefficient not found: ", coef_name,
      "\nAvailable coefficients:\n - ",
      paste(colnames(design), collapse = "\n - ")
    )
  }

  timing <- system.time({
    message("[DRIMSeq] ", tp, ": estimating precision")
    set.seed(123)
    d <- DRIMSeq::dmPrecision(d, design = design, BPPARAM = BPPARAM)

    message("[DRIMSeq] ", tp, ": fitting Dirichlet-multinomial model")
    d <- DRIMSeq::dmFit(d, design = design, BPPARAM = BPPARAM)

    message("[DRIMSeq] ", tp, ": performing likelihood ratio test")
    d <- DRIMSeq::dmTest(d, coef = coef_name)
  })

  message("[DRIMSeq] ", tp, ": runtime = ", round(timing["elapsed"], 2), " sec")

  res_gene <- as.data.frame(DRIMSeq::results(d))
  res_tx <- as.data.frame(DRIMSeq::results(d, level = "feature"))

  pScreen <- res_gene$pvalue
  pScreen[is.na(pScreen)] <- 1
  names(pScreen) <- res_gene$gene_id

  pConfirmation <- matrix(res_tx$pvalue, ncol = 1)
  pConfirmation[is.na(pConfirmation)] <- 1
  rownames(pConfirmation) <- res_tx$feature_id

  tx2g_stageR <- res_tx[, c("feature_id", "gene_id")]

  message("[DRIMSeq] ", tp, ": running stageR correction")

  st <- stageR::stageRTx(
    pScreen = pScreen,
    pConfirmation = pConfirmation,
    pScreenAdjusted = FALSE,
    tx2gene = tx2g_stageR
  )

  st <- stageR::stageWiseAdjustment(st, method = "dtu", alpha = padj_cutoff)

  gene_full <- .standardize_drimseq_gene_res(
    res_df = res_gene,
    timepoint = tp,
    contrast = cmp,
    annot_df = annot_df
  )

  tx_full <- .standardize_drimseq_tx_res(
    res_df = res_tx,
    timepoint = tp,
    contrast = cmp,
    annot_df = annot_df
  )

  sig_gene_ids <- suppressMessages(
    rownames(as.data.frame(stageR::getSignificantGenes(st)))
  )
  sig_tx_ids <- suppressMessages(
    rownames(as.data.frame(stageR::getSignificantTx(st)))
  )

  gene_full$stageR_sig <- gene_full$gene_id_full %in% sig_gene_ids
  gene_full$significant <- gene_full$stageR_sig

  tx_full$stageR_sig <- tx_full$feature_id %in% sig_tx_ids
  tx_full$significant <- tx_full$stageR_sig

  gene_sig <- gene_full[gene_full$stageR_sig, , drop = FALSE]
  gene_sig <- gene_sig[!duplicated(gene_sig$gene_id), , drop = FALSE]
  if ("padj" %in% colnames(gene_sig)) {
    gene_sig <- gene_sig[order(gene_sig$padj), , drop = FALSE]
  }

  tx_sig <- tx_full[tx_full$stageR_sig, , drop = FALSE]
  tx_sig <- tx_sig[!duplicated(tx_sig$feature_id), , drop = FALSE]
  if ("padj" %in% colnames(tx_sig)) {
    tx_sig <- tx_sig[order(tx_sig$padj), , drop = FALSE]
  }

  list(
    coef = coef_name,
    dm = d,
    stageR = st,
    counts = counts_df,
    results = list(
      gene_full = gene_full,
      gene_sig = gene_sig,
      tx_full = tx_full,
      tx_sig = tx_sig
    ),
    filtered_genes = unique(.strip_ens_version(DRIMSeq::counts(d)$gene_id))
  )
}

.import_drimseq_counts_all <- function(coldata, tx2gene) {
  files_named <- coldata$files
  names(files_named) <- coldata$names

  txi <- tximport::tximport(
    files_named,
    type = "salmon",
    tx2gene = tx2gene[, c("TXNAME", "GENEID")],
    txOut = TRUE,
    countsFromAbundance = "dtuScaledTPM"
  )

  mat <- txi$counts
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]

  map <- tx2gene[match(rownames(mat), tx2gene$TXNAME), c("GENEID", "TXNAME")]
  if (anyNA(map$GENEID) || anyNA(map$TXNAME)) {
    stop("Failed to map some transcripts back to tx2gene.")
  }

  data.frame(
    gene_id = map$GENEID,
    feature_id = map$TXNAME,
    mat,
    check.names = FALSE
  )
}

# -------------------------
# Public API
# -------------------------

run_drimseq <- function(coldata_tsv = NULL,
                        env_file = NULL,
                        resources_dir = NULL,
                        indexDir = NULL,
                        timepoints = c("H1", "H3", "H24"),
                        contrast_var = "condition",
                        ref_level = "C1",
                        padj_cutoff = 0.10,
                        annot_df = NULL,
                        filter_params = list(
                          C1 = list(
                            min_samps_feature_expr = 3, min_samps_gene_expr = 6,
                            min_feature_expr = 10, min_gene_expr = 10
                          ),
                          H1 = list(
                            min_samps_feature_expr = 3, min_samps_gene_expr = 6,
                            min_feature_expr = 10, min_gene_expr = 10
                          ),
                          H3 = list(
                            min_samps_feature_expr = 3, min_samps_gene_expr = 6,
                            min_feature_expr = 10, min_gene_expr = 10
                          ),
                          H24 = list(
                            min_samps_feature_expr = 2, min_samps_gene_expr = 5,
                            min_feature_expr = 10, min_gene_expr = 10
                          )
                        ),
                        drop_single_tx_genes = FALSE,
                        workers = NULL,
                        force = FALSE,
                        out_rds = NULL,
                        out_xlsx = NULL) {
  .check_pkg("tximport")
  .check_pkg("DRIMSeq")
  .check_pkg("stageR")
  .check_pkg("BiocParallel")
  .check_pkg("openxlsx")

  outdir <- file.path(.get_results_dir(), "analysis/drimseq")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(out_rds)) out_rds <- file.path(outdir, paste0("drimseq_results.", ref_level, ".rds"))
  if (is.null(out_xlsx)) out_xlsx <- file.path(outdir, paste0("drimseq_results.", ref_level, ".xlsx"))

  if (file.exists(out_rds) && !force) {
    message("Loading cached DRIMSeq results: ", out_rds)
    return(readRDS(out_rds))
  }

  if (is.null(coldata_tsv)) {
    coldata_tsv <- file.path(PROJECT_ROOT, "results/counts/salmon/coldata.tsv")
  }

  if (is.null(env_file)) {
    env_file <- file.path(PROJECT_ROOT, "config/environment/analysis_environment.env")
  }

  if (is.null(resources_dir)) {
    resources_dir <- file.path(PROJECT_ROOT, "resources")
  }

  if (is.null(indexDir)) {
    indexDir <- file.path(PROJECT_ROOT, "resources/salmon_index")
  }

  if (is.null(annot_df)) {
    annot_df <- .get_annot()
  }

  if (is.null(workers)) {
    workers <- .resolve_threads_from_env(env_file)
  }
  if (is.null(workers) || is.na(workers)) workers <- 1

  coldata <- .read_coldata_tsv(
    coldata_tsv = coldata_tsv,
    contrast_var = contrast_var,
    sample_col = "names",
    files_col = "files"
  )

  tx2gene <- .get_tx2gene(
    env_file = env_file,
    resources_dir = resources_dir,
    force = force
  )

  counts_all <- .import_drimseq_counts_all(
    coldata = coldata,
    tx2gene = tx2gene
  )

  BPPARAM <- if (.Platform$OS.type == "windows") {
    BiocParallel::SnowParam(workers = workers, progressbar = TRUE)
  } else {
    BiocParallel::MulticoreParam(workers = workers, progressbar = TRUE)
  }

  padj_pct <- formatC(padj_cutoff * 100, format = "f", digits = 0)
  message(
    "stageR notice: The returned adjusted p-values are based on a stage-wise testing ",
    "approach and are only valid for the provided target OFDR level of ",
    padj_pct, "%. If a different target OFDR level is of interest, ",
    "the entire adjustment should be re-run."
  )

  by_tp <- setNames(lapply(timepoints, function(tp) {
    if (!tp %in% names(filter_params)) {
      stop("No filter_params entry found for timepoint: ", tp)
    }

    .run_one_drimseq_timepoint(
      tp = tp,
      coldata = coldata,
      counts_all = counts_all,
      contrast_var = contrast_var,
      ref_level = ref_level,
      filter_params = filter_params[[tp]],
      padj_cutoff = padj_cutoff,
      drop_single_tx_genes = drop_single_tx_genes,
      BPPARAM = BPPARAM,
      annot_df = annot_df
    )
  }), timepoints)

  gene_full_list <- lapply(by_tp, function(x) x$results$gene_full)
  gene_sig_list <- lapply(by_tp, function(x) x$results$gene_sig)
  tx_full_list <- lapply(by_tp, function(x) x$results$tx_full)
  tx_sig_list <- lapply(by_tp, function(x) x$results$tx_sig)

  summary_tbl <- .make_drimseq_summary(
    gene_sig_list = gene_sig_list,
    tx_sig_list = tx_sig_list,
    padj_cutoff = padj_cutoff
  )

  .write_drimseq_xlsx(
    gene_sig_list = gene_sig_list,
    tx_sig_list = tx_sig_list,
    summary_tbl = summary_tbl,
    out_xlsx = out_xlsx
  )

  out <- list(
    meta = list(
      project_root = PROJECT_ROOT,
      coldata_tsv = coldata_tsv,
      env_file = env_file,
      resources_dir = resources_dir,
      indexDir = indexDir,
      contrast_var = contrast_var,
      ref_level = ref_level,
      timepoints = timepoints,
      padj_cutoff = padj_cutoff,
      filter_params = filter_params,
      workers = workers
    ),
    common_data = list(
      coldata = coldata,
      counts_all = counts_all,
      tx2gene = tx2gene
    ),
    by_tp = by_tp,
    results = list(
      gene_full_by_tp = gene_full_list,
      gene_sig_by_tp = gene_sig_list,
      tx_full_by_tp = tx_full_list,
      tx_sig_by_tp = tx_sig_list,
      gene_full_all = do.call(rbind, gene_full_list),
      gene_sig_all = do.call(rbind, gene_sig_list),
      tx_full_all = do.call(rbind, tx_full_list),
      tx_sig_all = do.call(rbind, tx_sig_list),
      sig_all = list(
        genes = do.call(rbind, gene_sig_list),
        transcripts = do.call(rbind, tx_sig_list)
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
