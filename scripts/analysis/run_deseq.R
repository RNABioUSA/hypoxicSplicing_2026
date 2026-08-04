PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
if (!nzchar(PROJECT_ROOT)) stop("PROJECT_ROOT env var not set.")

HELPERS_FILE <- file.path(PROJECT_ROOT, "scripts/utils/helpers.R")
if (!file.exists(HELPERS_FILE)) stop("helpers.R not found at: ", HELPERS_FILE)
source(HELPERS_FILE)

# ------------------------------------------------------------------------------
# Internal Helpers
# ------------------------------------------------------------------------------
.standardize_deseq2_res <- function(res_df,
                                    timepoint,
                                    contrast,
                                    tool = "deseq2",
                                    annot_df = NULL,
                                    padj_cutoff = 0.10) {
  out <- res_df
  out$ensgene <- rownames(out)
  rownames(out) <- NULL

  out$ensgene <- .strip_ens_version(out$ensgene)

  out$timepoint <- timepoint
  out$contrast <- contrast
  out$tool <- tool
  out$comparison <- contrast

  out <- .attach_gene_annot(out, annot_df = annot_df)

  out$direction <- ifelse(
    is.na(out$padj),
    NA_character_,
    ifelse(out$padj < padj_cutoff & out$log2FoldChange > 0, "up",
      ifelse(out$padj < padj_cutoff & out$log2FoldChange < 0, "down", "ns")
    )
  )

  out$neglog10_padj <- ifelse(
    is.na(out$padj) | out$padj <= 0,
    NA_real_,
    -log10(out$padj)
  )

  out
}

.run_one_timepoint <- function(dds,
                               timepoint,
                               ref_level,
                               contrast_var = "condition",
                               shrink_type = "ashr",
                               padj_cutoff = 0.10,
                               annot_df = NULL) {
  coef_name <- sprintf("%s_%s_vs_%s", contrast_var, timepoint, ref_level)

  rn <- DESeq2::resultsNames(dds)
  if (!(coef_name %in% rn)) {
    stop(
      "Coefficient not found: ", coef_name,
      "\nAvailable:\n  - ", paste(rn, collapse = "\n  - ")
    )
  }

  res <- DESeq2::lfcShrink(dds, coef = coef_name, type = shrink_type)

  full <- .standardize_deseq2_res(
    res_df = as.data.frame(res),
    timepoint = timepoint,
    contrast = coef_name,
    annot_df = annot_df,
    padj_cutoff = padj_cutoff
  )

  sig <- full[!duplicated(full$ensgene) & !is.na(full$padj) & full$padj < padj_cutoff, , drop = FALSE]
  sig <- sig[order(sig$padj), , drop = FALSE]

  list(coef = coef_name, full = full, sig = sig)
}

.make_summary <- function(sig_list, padj_cutoff = 0.10) {
  do.call(rbind, lapply(names(sig_list), function(tp) {
    df <- sig_list[[tp]]
    if (nrow(df) == 0) {
      data.frame(
        timepoint = tp,
        padj_cutoff = padj_cutoff,
        n_sig = 0L,
        n_up = 0L,
        n_down = 0L
      )
    } else {
      data.frame(
        timepoint = tp,
        padj_cutoff = padj_cutoff,
        n_sig = nrow(df),
        n_up = sum(df$direction == "up", na.rm = TRUE),
        n_down = sum(df$direction == "down", na.rm = TRUE)
      )
    }
  }))
}

.write_sig_xlsx <- function(sig_list, summary_tbl, out_xlsx) {
  .check_pkg("openxlsx")

  wb <- openxlsx::createWorkbook()

  openxlsx::addWorksheet(wb, "Summary")
  openxlsx::writeData(wb, "Summary", summary_tbl)

  for (tp in names(sig_list)) {
    openxlsx::addWorksheet(wb, tp)
    df <- sig_list[[tp]]

    if (nrow(df) == 0) {
      openxlsx::writeData(wb, tp, data.frame(note = "No significant genes at this cutoff."))
    } else {
      openxlsx::writeData(wb, tp, df)
    }
  }

  dir.create(dirname(out_xlsx), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
  invisible(TRUE)
}

.get_deseq2_tested_universe <- function(full_list) {
  lapply(full_list, function(df) {
    unique(df$ensgene[!is.na(df$padj)])
  })
}

# ------------------------------------------------------------------------------
# Public API
# ------------------------------------------------------------------------------
run_deseq <- function(coldata_tsv = NULL,
                      env_file = NULL,
                      resources_dir = NULL,
                      indexDir = NULL,
                      timepoints = c("H1", "H3", "H24"),
                      design = ~condition,
                      contrast_var = "condition",
                      ref_level = "C1",
                      padj_cutoff = 0.10,
                      shrink_type = "ashr",
                      assignRanges = "abundant",
                      annot_df = NULL,
                      force = FALSE,
                      out_rds = NULL,
                      out_xlsx = NULL) {
  .check_pkg("tximeta")
  .check_pkg("DESeq2")
  .check_pkg("ashr")

  outdir <- file.path(.get_results_dir(), "analysis", "deseq")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(out_rds)) {
    out_rds <- file.path(outdir, paste0("deseq_results.", ref_level, ".rds"))
  }
  if (is.null(out_xlsx)) {
    out_xlsx <- file.path(outdir, paste0("deseq_results.", ref_level, ".xlsx"))
  }

  if (file.exists(out_rds) && !force) {
    message("Loading cached DESeq2 results: ", out_rds)
    return(readRDS(out_rds))
  }

  paths <- .get_runtime_paths()

  coldata_tsv <- coldata_tsv %||% paths$coldata_tsv
  env_file <- env_file %||% paths$env_file
  resources_dir <- resources_dir %||% paths$resources_dir
  indexDir <- indexDir %||% paths$indexDir

  if (is.null(annot_df)) {
    annot_df <- .get_annot()
  }

  coldata <- .read_coldata_tsv(
    coldata_tsv = coldata_tsv,
    contrast_var = contrast_var,
    sample_col = "names",
    files_col = "files"
  )

  coldata[[contrast_var]] <- factor(coldata[[contrast_var]])
  if (!(ref_level %in% levels(coldata[[contrast_var]]))) {
    stop("Reference level not found in ", contrast_var, ": ", ref_level)
  }
  coldata[[contrast_var]] <- stats::relevel(coldata[[contrast_var]], ref = ref_level)

  gse <- .get_gse(
    coldata_tsv = coldata_tsv,
    env_file = env_file,
    resources_dir = resources_dir,
    indexDir = indexDir,
    force = force,
    assignRanges = assignRanges
  )

  dds <- DESeq2::DESeqDataSet(gse, design = design)

  SummarizedExperiment::colData(dds)[[contrast_var]] <- stats::relevel(
    factor(SummarizedExperiment::colData(dds)[[contrast_var]]),
    ref = ref_level
  )

  dds <- DESeq2::DESeq(dds)

  by_tp <- setNames(lapply(timepoints, function(tp) {
    .run_one_timepoint(
      dds = dds,
      timepoint = tp,
      ref_level = ref_level,
      contrast_var = contrast_var,
      shrink_type = shrink_type,
      padj_cutoff = padj_cutoff,
      annot_df = annot_df
    )
  }), timepoints)

  by_tp <- lapply(names(by_tp), function(tp) {
    obj <- by_tp[[tp]]
    cmp <- paste0(tp, "_vs_", ref_level)

    obj$full$comparison <- cmp
    obj$full$gene_id <- obj$full$ensgene

    obj$sig$comparison <- cmp
    obj$sig$gene_id <- obj$sig$ensgene

    obj
  }) |>
    setNames(timepoints)

  full_list <- lapply(by_tp, `[[`, "full")
  sig_list <- lapply(by_tp, `[[`, "sig")
  tested_universe_by_tp <- .get_deseq2_tested_universe(full_list)
  tested_universe_all <- unique(unlist(tested_universe_by_tp))

  summary_tbl <- .make_summary(sig_list, padj_cutoff = padj_cutoff)
  .write_sig_xlsx(sig_list, summary_tbl, out_xlsx)

  out <- list(
    meta = list(
      project_root = PROJECT_ROOT,
      env_file = env_file,
      resources_dir = resources_dir,
      indexDir = indexDir,
      coldata_tsv = coldata_tsv,
      design = deparse(design),
      contrast_var = contrast_var,
      ref_level = ref_level,
      timepoints = timepoints,
      padj_cutoff = padj_cutoff,
      shrink_type = shrink_type,
      assignRanges = assignRanges
    ),
    dds = dds,
    by_tp = by_tp,
    results = list(
      gene_full_by_tp = full_list,
      gene_sig_by_tp = sig_list,
      gene_full_all = do.call(rbind, full_list),
      gene_sig_all = do.call(rbind, sig_list),
      tested_universe_by_tp = tested_universe_by_tp,
      tested_universe_all = tested_universe_all
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
