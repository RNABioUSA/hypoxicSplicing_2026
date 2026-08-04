# ============================================================
# plot_sup_fig_4.R
# ============================================================
# Purpose:
#   Evaluate robustness of the observed limited overlap between
#   gene-level differential expression and splicing regulation.
#
#   Specifically tests whether low concordance is stable across:
#     1. gene-level statistical method (DESeq2, edgeR, limma-voom)
#     2. significance threshold (FDR)
#
#   All analyses are performed on a harmonized input universe
#   defined by DRIMSeq filtering, ensuring that gene- and
#   splicing-level methods are evaluated on matched input sets.
#
# Harmonization Strategy:
#   - DRIMSeq defines the gene and transcript universe per timepoint.
#   - Gene-level methods (DESeq2, edgeR, limma-voom) are rerun on
#     DRIMSeq-filtered count data for each pairwise comparison.
#   - Splicing methods:
#       * DEXSeq is rerun on filtered SummarizedExperiment objects.
#       * SUPPA2 is rerun using filtered isoform TPM inputs.
#       * DRIMSeq results inherently reflect the filtered universe.
#
#   This ensures that all overlap comparisons reflect true biological
#   differences rather than differences in input filtering.
#
# Notes:
#   - No post hoc restriction of standard results is used in the
#     primary analysis; all results are derived from harmonized inputs.
#   - Threshold sensitivity is evaluated using harmonized results only.
#   - All gene sets are validated against the DRIMSeq-derived universe
#     prior to plotting.

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
if (!nzchar(PROJECT_ROOT)) stop("PROJECT_ROOT env var not set.")

HELPERS_FILE <- file.path(PROJECT_ROOT, "scripts/utils/helpers.R")
if (!file.exists(HELPERS_FILE)) stop("helpers.R not found at: ", HELPERS_FILE)
source(HELPERS_FILE)

COLORS_FILE <- file.path(PROJECT_ROOT, "scripts/utils/color_palette.R")
if (!file.exists(COLORS_FILE)) stop("color_palette.R not found at: ", COLORS_FILE)
source(COLORS_FILE)

# --------------------------------------------------
# Internal Helpers
# --------------------------------------------------

.timepoint_order <- function() {
  c("H1", "H3", "H24")
}

.timepoint_labels <- function() {
  stats::setNames(
    names(timepoint_base_colors)[grepl("Hypoxia", names(timepoint_base_colors))],
    .timepoint_order()
  )
}

.timepoint_shades <- function(timepoint) {
  .get_timepoint_shaded_colors(
    timepoints = timepoint,
    labels = c("Transcript Level", "Overlap", "Gene Level")
  )[[timepoint]]
}

.timepoint_main_color <- function(timepoint) {
  .timepoint_shades(timepoint)[["Overlap"]]
}

.timepoint_layer_colors <- function(timepoint) {
  .timepoint_shades(timepoint)[c("Gene Level", "Transcript Level", "Overlap")]
}

.overlap_category_order <- function() {
  c("Transcript Level", "Overlap", "Gene Level")
}

.tool_order <- function() {
  c(
    .splicing_methods(),
    .gene_methods()
  )
}

.tool_colors <- function(shade = "base") {
  base_cols <- c(
    "DESeq2" = unname(category_base_colors[["Category 2"]]),
    "edgeR" = unname(category_base_colors[["Category 4"]]),
    "limma-voom" = unname(category_base_colors[["Category 5"]]),
    "DEXSeq" = unname(splicing_tool_colors[["DEXSeq"]]),
    "DRIMSeq" = unname(splicing_tool_colors[["DRIMSeq"]]),
    "SUPPA2" = unname(splicing_tool_colors[["SUPPA2"]])
  )

  shaded <- .generate_shaded_palette(
    base_colors = base_cols,
    labels = c("light", "base", "dark")
  )

  stats::setNames(
    vapply(shaded, function(x) unname(x[[shade]]), character(1)),
    names(shaded)
  )
}

.gene_methods <- function() {
  c("DESeq2", "edgeR", "limma-voom")
}

.splicing_methods <- function() {
  c("DEXSeq", "DRIMSeq", "SUPPA2")
}

.get_timepoint_min_size <- function(min_size, timepoint) {
  if (length(min_size) == 1 && is.null(names(min_size))) {
    return(min_size)
  }

  if (is.null(names(min_size))) {
    stop("If min_size has length > 1, it must be named by timepoint.")
  }

  if (!timepoint %in% names(min_size)) {
    stop("Missing min_size for timepoint: ", timepoint)
  }

  min_size[[timepoint]]
}

.require_result_columns <- function(df, cols, name = deparse(substitute(df))) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(name, " missing required column(s): ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

.as_gene_set_tbl <- function(genes, timepoint, method, layer, source = NA_character_) {
  genes <- unique(.strip_ens_version(genes))
  genes <- genes[!is.na(genes) & genes != ""]

  data.frame(
    timepoint = timepoint,
    method = method,
    layer = layer,
    source = source,
    ensgene = genes,
    stringsAsFactors = FALSE
  )
}

# --------------------------------------------------
# Harmonized Universe Helpers
# --------------------------------------------------

.get_harmonized_universe_from_drimseq <- function(drimseq_results) {
  .check_pkg("dplyr")

  if (is.null(drimseq_results$by_tp)) {
    stop("Expected drimseq_results$by_tp.")
  }

  out <- lapply(.timepoint_order(), function(tp) {
    x <- drimseq_results$by_tp[[tp]]
    if (is.null(x)) stop("DRIMSeq by_tp missing timepoint: ", tp)

    # Current DRIMSeq scripts store filtered genes directly.
    gene_ids <- x$filtered_genes

    # Fallback: derive from retained DRIMSeq counts if available.
    if (is.null(gene_ids) && !is.null(x$counts$gene_id)) {
      gene_ids <- unique(x$counts$gene_id)
    }

    if (is.null(gene_ids)) {
      stop(
        "Could not derive harmonized gene universe for ", tp,
        ". Expected by_tp[[tp]]$filtered_genes or by_tp[[tp]]$counts$gene_id."
      )
    }

    data.frame(
      timepoint = tp,
      ensgene = unique(.strip_ens_version(gene_ids)),
      stringsAsFactors = FALSE
    )
  }) |>
    dplyr::bind_rows() |>
    dplyr::distinct(timepoint, ensgene)

  invisible(out)
}

.get_harmonized_tx_universe_from_drimseq <- function(drimseq_results) {
  .check_pkg("dplyr")

  if (is.null(drimseq_results$by_tp)) {
    stop("Expected drimseq_results$by_tp.")
  }

  out <- lapply(.timepoint_order(), function(tp) {
    x <- drimseq_results$by_tp[[tp]]
    if (is.null(x)) stop("DRIMSeq by_tp missing timepoint: ", tp)

    tx_ids <- NULL

    if (!is.null(x$filtered_tx)) {
      tx_ids <- x$filtered_tx
    } else if (!is.null(x$counts$feature_id)) {
      tx_ids <- unique(x$counts$feature_id)
    } else if (!is.null(x$dmData)) {
      tx_ids <- unique(DRIMSeq::counts(x$dmData)$feature_id)
    }

    if (is.null(tx_ids)) {
      stop(
        "Could not derive transcript universe for ", tp,
        ". Expected filtered_tx, counts$feature_id, or dmData."
      )
    }

    data.frame(
      timepoint = tp,
      enstx = unique(.strip_ens_version(tx_ids)),
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(out) |>
    dplyr::distinct(timepoint, enstx)
}

.subset_gse_for_pairwise_universe <- function(gse, timepoint, universe_genes, ref_level = "C1") {
  keep_samples <- SummarizedExperiment::colData(gse)$condition %in% c(ref_level, timepoint)
  gse_pair <- gse[, keep_samples]

  gene_ids <- .strip_ens_version(rownames(gse_pair))
  keep_genes <- gene_ids %in% universe_genes
  gse_pair[keep_genes, , drop = FALSE]
}

.run_deseq2_pairwise_harmonized <- function(
  gse_obj,
  timepoint,
  ref_level = "C1",
  fdr = 0.10,
  shrink_type = "apeglm"
) {
  .check_pkg(c("DESeq2", "dplyr", "tibble"))

  dds <- suppressWarnings(DESeq2::DESeqDataSet(gse_obj, design = ~condition))
  dds$condition <- stats::relevel(factor(dds$condition), ref = ref_level)
  dds <- DESeq2::DESeq(dds, minReplicatesForReplace = Inf)

  res <- DESeq2::results(
    dds,
    contrast = c("condition", timepoint, ref_level),
    cooksCutoff = FALSE,
    independentFiltering = FALSE
  )

  rn <- DESeq2::resultsNames(dds)
  coef_name <- rn[grepl(paste0("^condition_", timepoint, "_vs_", ref_level, "$"), rn)]
  if (length(coef_name) != 1) {
    coef_name <- rn[grepl("^condition_", rn)]
  }
  if (length(coef_name) != 1) {
    stop("Could not uniquely identify DESeq2 coefficient. Available: ", paste(rn, collapse = ", "))
  }

  res_shr <- DESeq2::lfcShrink(dds, coef = coef_name, type = shrink_type, res = res)

  full <- as.data.frame(res_shr) |>
    tibble::rownames_to_column("ensgene") |>
    dplyr::mutate(
      ensgene = .strip_ens_version(ensgene),
      timepoint = timepoint,
      method = "DESeq2",
      layer = "Gene",
      significant = !is.na(padj) & padj < fdr
    )

  list(
    full = full,
    sig = dplyr::filter(full, significant),
    dds = dds,
    coef = coef_name
  )
}

.run_edger_pairwise_harmonized <- function(
  gse_obj,
  timepoint,
  ref_level = "C1",
  fdr = 0.10
) {
  .check_pkg(c("edgeR", "limma", "dplyr", "tibble"))

  y <- edgeR::DGEList(counts = SummarizedExperiment::assays(gse_obj)[["counts"]])
  grp <- factor(SummarizedExperiment::colData(gse_obj)$condition, levels = c(ref_level, timepoint))
  design <- stats::model.matrix(~ 0 + grp)
  colnames(design) <- levels(grp)

  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmQLFit(y, design)
  contr <- limma::makeContrasts(contrasts = paste0(timepoint, "-", ref_level), levels = design)
  qlf <- edgeR::glmQLFTest(fit, contrast = contr)

  full <- as.data.frame(edgeR::topTags(qlf, n = Inf)$table) |>
    tibble::rownames_to_column("ensgene") |>
    dplyr::mutate(
      ensgene = .strip_ens_version(ensgene),
      timepoint = timepoint,
      method = "edgeR",
      layer = "Gene",
      padj = FDR,
      significant = !is.na(FDR) & FDR < fdr
    )

  list(
    full = full,
    sig = dplyr::filter(full, significant),
    fit = fit,
    contrast = contr
  )
}

.run_voom_pairwise_harmonized <- function(
  gse_obj,
  timepoint,
  ref_level = "C1",
  fdr = 0.10
) {
  .check_pkg(c("edgeR", "limma", "dplyr", "tibble"))

  y <- edgeR::DGEList(counts = SummarizedExperiment::assays(gse_obj)[["counts"]])
  grp <- factor(SummarizedExperiment::colData(gse_obj)$condition, levels = c(ref_level, timepoint))
  design <- stats::model.matrix(~ 0 + grp)
  colnames(design) <- levels(grp)

  y <- edgeR::calcNormFactors(y)
  v <- limma::voom(y, design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  contr <- limma::makeContrasts(contrasts = paste0(timepoint, "-", ref_level), levels = design)
  fit2 <- limma::contrasts.fit(fit, contr)
  fit2 <- limma::eBayes(fit2)

  full <- as.data.frame(limma::topTable(fit2, n = Inf, adjust.method = "BH")) |>
    tibble::rownames_to_column("ensgene") |>
    dplyr::mutate(
      ensgene = .strip_ens_version(ensgene),
      timepoint = timepoint,
      method = "limma-voom",
      layer = "Gene",
      padj = adj.P.Val,
      significant = !is.na(adj.P.Val) & adj.P.Val < fdr
    )

  list(
    full = full,
    sig = dplyr::filter(full, significant),
    fit = fit2
  )
}

.prepare_suppa_harmonized_isotpm <- function(
  drimseq_results,
  isotpm_dir = file.path(.get_results_dir(), "analysis", "suppa", "standard", "C1", "iso_tpm"),
  outdir = file.path(.get_results_dir(), "analysis", "suppa", "harmonized", "filtered_inputs")
) {
  tx_universe <- .get_harmonized_tx_universe_from_drimseq(drimseq_results)

  iso_paths <- list(
    C1  = file.path(isotpm_dir, "iso_tpm_C1.txt"),
    H1  = file.path(isotpm_dir, "iso_tpm_H1.txt"),
    H3  = file.path(isotpm_dir, "iso_tpm_H3.txt"),
    H24 = file.path(isotpm_dir, "iso_tpm_H24.txt")
  )

  iso <- lapply(iso_paths, .read_suppa_isotpm)

  out <- lapply(.timepoint_order(), function(tp) {
    tx_tp <- tx_universe$enstx[tx_universe$timepoint == tp]

    ref_filt <- .filter_suppa_isotpm_by_tx(iso$C1$df, tx_tp)
    test_filt <- .filter_suppa_isotpm_by_tx(iso[[tp]]$df, tx_tp)

    ref_path <- file.path(outdir, paste0("iso_tpm_C1v", tp, "_filtered.txt"))
    test_path <- file.path(outdir, paste0("iso_tpm_", tp, "_filtered.txt"))

    .write_suppa_isotpm(
      ref_filt,
      ref_path,
      sample_names = iso$C1$sample_names
    )

    .write_suppa_isotpm(
      test_filt,
      test_path,
      sample_names = iso[[tp]]$sample_names
    )

    data.frame(
      timepoint = tp,
      ref_path = ref_path,
      test_path = test_path,
      n_ref_tx = nrow(ref_filt),
      n_test_tx = nrow(test_filt),
      stringsAsFactors = FALSE
    )
  })

  list(
    outdir = outdir,
    tx_universe = tx_universe,
    isotpm_paths = dplyr::bind_rows(out)
  )
}

.extract_splicing_sets_from_results <- function(
  dexseq_results = NULL,
  drimseq_results = NULL,
  suppa_results = NULL,
  thresholds = c(0.10, 0.05, 0.01),
  universe = NULL,
  restrict_to_universe = TRUE,
  source = "harmonized_splicing_input"
) {
  .check_pkg("dplyr")

  out <- list()

  if (!is.null(dexseq_results)) {
    dex_list <- dexseq_results$results$gene_full_by_tp
    if (is.null(dex_list)) {
      stop("DEXSeq object must contain results$gene_full_by_tp.")
    }

    out$dexseq <- dplyr::bind_rows(lapply(names(dex_list), function(tp) {
      df <- dex_list[[tp]]
      .require_result_columns(df, c("ensgene", "padj"), paste0("DEXSeq ", tp))

      dplyr::bind_rows(lapply(thresholds, function(thr) {
        genes <- unique(.strip_ens_version(
          df$ensgene[!is.na(df$padj) & df$padj < thr]
        ))

        if (isTRUE(restrict_to_universe) && !is.null(universe)) {
          genes <- intersect(genes, universe$ensgene[universe$timepoint == tp])
        }

        .as_gene_set_tbl(genes, tp, "DEXSeq", "Splicing", source) |>
          dplyr::mutate(fdr = thr)
      }))
    }))
  }

  if (!is.null(drimseq_results)) {
    drim_list <- drimseq_results$results$gene_full_by_tp
    if (is.null(drim_list)) {
      stop("DRIMSeq object must contain results$gene_full_by_tp.")
    }

    out$drimseq <- dplyr::bind_rows(lapply(names(drim_list), function(tp) {
      df <- drim_list[[tp]]
      .require_result_columns(df, c("ensgene", "padj"), paste0("DRIMSeq ", tp))

      dplyr::bind_rows(lapply(thresholds, function(thr) {
        genes <- unique(.strip_ens_version(
          df$ensgene[!is.na(df$padj) & df$padj < thr]
        ))

        if (isTRUE(restrict_to_universe) && !is.null(universe)) {
          genes <- intersect(genes, universe$ensgene[universe$timepoint == tp])
        }

        .as_gene_set_tbl(genes, tp, "DRIMSeq", "Splicing", source) |>
          dplyr::mutate(fdr = thr)
      }))
    }))
  }

  if (!is.null(suppa_results)) {
    suppa_list <- suppa_results$results$event_full_by_tp
    if (is.null(suppa_list)) {
      stop("SUPPA2 object must contain results$event_full_by_tp.")
    }

    out$suppa2 <- dplyr::bind_rows(lapply(names(suppa_list), function(tp) {
      df <- suppa_list[[tp]]
      .require_result_columns(df, c("ensgene", "pvalue"), paste0("SUPPA2 ", tp))

      dplyr::bind_rows(lapply(thresholds, function(thr) {
        genes <- unique(.strip_ens_version(
          df$ensgene[!is.na(df$pvalue) & df$pvalue < thr]
        ))

        if (isTRUE(restrict_to_universe) && !is.null(universe)) {
          genes <- intersect(genes, universe$ensgene[universe$timepoint == tp])
        }

        .as_gene_set_tbl(genes, tp, "SUPPA2", "Splicing", source) |>
          dplyr::mutate(fdr = thr)
      }))
    }))
  }

  dplyr::bind_rows(out)
}

.extract_gene_sets_from_harmonized_results <- function(
  harmonized_gene_results,
  thresholds = c(0.10, 0.05, 0.01),
  methods = .gene_methods()
) {
  .check_pkg("dplyr")

  p_col_by_method <- c(
    "DESeq2" = "padj",
    "edgeR" = "FDR",
    "limma-voom" = "adj.P.Val"
  )

  dplyr::bind_rows(lapply(names(harmonized_gene_results$by_tp), function(tp) {
    dplyr::bind_rows(lapply(methods, function(method) {
      full <- harmonized_gene_results$by_tp[[tp]][[method]]$full
      p_col <- p_col_by_method[[method]]

      .require_result_columns(full, c("ensgene", p_col), paste0(method, " ", tp))

      dplyr::bind_rows(lapply(thresholds, function(thr) {
        genes <- unique(.strip_ens_version(
          full$ensgene[!is.na(full[[p_col]]) & full[[p_col]] < thr]
        ))

        .as_gene_set_tbl(
          genes = genes,
          timepoint = tp,
          method = method,
          layer = "Gene",
          source = paste0("harmonized_gene_input_fdr_", thr)
        ) |>
          dplyr::mutate(fdr = thr)
      }))
    }))
  }))
}

# --------------------------------------------------
# Main Analysis Object Builders
# --------------------------------------------------

make_filter_threshold_robustness_data <- function(
  dexseq_results,
  drimseq_results,
  suppa_results,
  harmonized_gene_results = NULL,
  thresholds = c(0.10, 0.05, 0.01),
  run_harmonized_gene_methods = TRUE,
  restrict_splicing_to_harmonized_universe = TRUE,
  ...
) {
  .check_pkg("dplyr")

  universe <- .get_harmonized_universe_from_drimseq(drimseq_results)

  if (is.null(harmonized_gene_results) && isTRUE(run_harmonized_gene_methods)) {
    harmonized_gene_results <- run_gene_methods_harmonized(
      drimseq_results = drimseq_results,
      fdr = max(thresholds),
      ...
    )
  }

  gene_sets_harmonized <- NULL
  if (!is.null(harmonized_gene_results)) {
    gene_sets_harmonized <- harmonized_gene_results$sig_tbl |>
      dplyr::mutate(fdr = max(thresholds))
  }

  # Threshold panel uses harmonized gene-level results.
  gene_sets_threshold <- .extract_gene_sets_from_harmonized_results(
    harmonized_gene_results = harmonized_gene_results,
    thresholds = thresholds,
    methods = .gene_methods()
  )

  splice_sets_threshold <- .extract_splicing_sets_from_results(
    dexseq_results = dexseq_results,
    drimseq_results = drimseq_results,
    suppa_results = suppa_results,
    thresholds = thresholds,
    universe = universe,
    restrict_to_universe = restrict_splicing_to_harmonized_universe,
    source = ifelse(
      restrict_splicing_to_harmonized_universe,
      "harmonized_splicing_input_restricted_to_harmonized_universe",
      "harmonized_splicing_input"
    )
  )

  # Harmonized-method UpSet panel at the primary FDR threshold.
  splice_sets_primary <- splice_sets_threshold |>
    dplyr::filter(fdr == max(thresholds))

  if (!is.null(gene_sets_harmonized)) {
    upset_sets <- dplyr::bind_rows(
      gene_sets_harmonized,
      splice_sets_primary
    )
  } else {
    upset_sets <- dplyr::bind_rows(
      gene_sets_threshold |> dplyr::filter(fdr == max(thresholds)),
      splice_sets_primary
    )
  }

  threshold_overlap <- compute_threshold_overlap_summary(
    gene_sets = gene_sets_threshold,
    splice_sets = splice_sets_threshold
  )

  out <- list(
    meta = list(
      thresholds = thresholds,
      restrict_splicing_to_harmonized_universe = restrict_splicing_to_harmonized_universe
    ),
    universe = universe,
    harmonized_gene_results = harmonized_gene_results,
    upset_sets = upset_sets,
    threshold_sets = list(
      gene = gene_sets_threshold,
      splicing = splice_sets_threshold
    ),
    threshold_overlap = threshold_overlap
  )

  out$harmonization_qc <- check_harmonized_analysis_universe(out)

  invisible(out)
}

compute_threshold_overlap_summary <- function(gene_sets, splice_sets) {
  .check_pkg("dplyr")

  fdr_levels <- sort(unique(c(gene_sets$fdr, splice_sets$fdr)), decreasing = FALSE)
  tps <- .timepoint_order()

  dplyr::bind_rows(lapply(tps, function(tp) {
    dplyr::bind_rows(lapply(fdr_levels, function(thr) {
      genes <- unique(gene_sets$ensgene[gene_sets$timepoint == tp & gene_sets$fdr == thr])
      splice <- unique(splice_sets$ensgene[splice_sets$timepoint == tp & splice_sets$fdr == thr])

      data.frame(
        timepoint = tp,
        fdr = thr,
        class = c("Splicing only", "Overlap", "Gene only"),
        n_genes = c(
          length(setdiff(splice, genes)),
          length(intersect(splice, genes)),
          length(setdiff(genes, splice))
        ),
        stringsAsFactors = FALSE
      )
    }))
  })) |>
    dplyr::group_by(timepoint, fdr) |>
    dplyr::mutate(prop = n_genes / sum(n_genes)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      timepoint = factor(timepoint, levels = .timepoint_order()),
      class = factor(class, levels = c("Splicing only", "Overlap", "Gene only")),
      fdr_label = factor(paste0("FDR < ", fdr, "\n"), levels = paste0("FDR < ", fdr_levels, "\n"))
    )
}

run_gene_methods_harmonized <- function(
  drimseq_results,
  coldata_tsv = NULL,
  env_file = NULL,
  resources_dir = NULL,
  indexDir = NULL,
  assignRanges = "abundant",
  ref_level = "C1",
  fdr = 0.10,
  force_gse = FALSE
) {
  .check_pkg(c("dplyr", "SummarizedExperiment"))

  paths <- .get_runtime_paths()
  if (is.null(coldata_tsv)) coldata_tsv <- paths$coldata_tsv
  if (is.null(env_file)) env_file <- paths$env_file
  if (is.null(resources_dir)) resources_dir <- paths$resources_dir
  if (is.null(indexDir)) indexDir <- paths$indexDir

  universe <- .get_harmonized_universe_from_drimseq(drimseq_results)

  gse <- .get_gse(
    coldata_tsv = coldata_tsv,
    env_file = env_file,
    resources_dir = resources_dir,
    indexDir = indexDir,
    force = force_gse,
    assignRanges = assignRanges
  )

  by_tp <- lapply(.timepoint_order(), function(tp) {
    univ_tp <- universe$ensgene[universe$timepoint == tp]
    gse_tp <- .subset_gse_for_pairwise_universe(gse, tp, univ_tp, ref_level = ref_level)

    list(
      DESeq2 = .run_deseq2_pairwise_harmonized(gse_tp, tp, ref_level = ref_level, fdr = fdr),
      edgeR = .run_edger_pairwise_harmonized(gse_tp, tp, ref_level = ref_level, fdr = fdr),
      `limma-voom` = .run_voom_pairwise_harmonized(gse_tp, tp, ref_level = ref_level, fdr = fdr),
      universe_n = length(unique(univ_tp))
    )
  })
  names(by_tp) <- .timepoint_order()

  sig_tbl <- dplyr::bind_rows(lapply(names(by_tp), function(tp) {
    dplyr::bind_rows(lapply(.gene_methods(), function(method) {
      .as_gene_set_tbl(
        genes = by_tp[[tp]][[method]]$sig$ensgene,
        timepoint = tp,
        method = method,
        layer = "Gene",
        source = "harmonized_gene_input"
      )
    }))
  }))

  invisible(list(
    meta = list(ref_level = ref_level, fdr = fdr),
    universe = universe,
    by_tp = by_tp,
    sig_tbl = sig_tbl
  ))
}

check_harmonized_analysis_universe <- function(plot_data, stop_on_fail = TRUE) {
  .check_pkg("dplyr")

  universe <- plot_data$universe
  sets <- dplyr::bind_rows(
    plot_data$threshold_sets$gene,
    plot_data$threshold_sets$splicing
  )

  qc <- sets |>
    dplyr::group_by(timepoint, layer, method, fdr) |>
    dplyr::summarise(
      n_detected = dplyr::n_distinct(ensgene),
      n_outside_universe = length(setdiff(
        unique(ensgene),
        universe$ensgene[universe$timepoint == dplyr::first(timepoint)]
      )),
      .groups = "drop"
    ) |>
    dplyr::mutate(pass = n_outside_universe == 0)

  if (any(!qc$pass)) {
    bad <- qc |> dplyr::filter(!pass)

    msg <- paste0(
      "Detected genes outside harmonized universe:\n",
      paste(
        bad$timepoint,
        bad$layer,
        bad$method,
        bad$fdr,
        "n_outside =",
        bad$n_outside_universe,
        collapse = "\n"
      )
    )

    if (isTRUE(stop_on_fail)) stop(msg)
    warning(msg)
  } else {
    message("[Sup Fig 4 QC] PASS: TRUE")
  }

  qc
}

# --------------------------------------------------
# Plotting Helpers
# --------------------------------------------------

.make_membership_matrix <- function(gene_sets, methods) {
  all_genes <- sort(unique(gene_sets$ensgene))
  mat <- data.frame(ensgene = all_genes, stringsAsFactors = FALSE)

  for (m in methods) {
    mat[[m]] <- mat$ensgene %in% gene_sets$ensgene[gene_sets$method == m]
  }

  mat
}

plot_harmonized_method_upset <- function(
  plot_data,
  timepoint,
  methods = .tool_order(),
  min_size = 10,
  base_size = 12
) {
  .check_pkg(c("dplyr", "ggplot2", "ComplexUpset"))

  df <- plot_data$upset_sets |>
    dplyr::filter(timepoint == !!timepoint)

  methods <- intersect(.tool_order(), unique(df$method))
  mat <- .make_membership_matrix(df, methods)

  tp_cols <- .timepoint_shades(timepoint)
  bar_col <- tp_cols[["Overlap"]]
  inactive_dot_col <- tp_cols[["Transcript Level"]]

  set_metadata <- data.frame(
    set = rev(methods),
    layer = ifelse(rev(methods) %in% .gene_methods(), "Gene Level", "Transcript Level")
  )

  set_metadata$layer <- factor(
    set_metadata$layer,
    levels = c("Gene Level", "Transcript Level")
  )

  stripe_cols <- tp_cols[c("Gene Level", "Transcript Level")]

  suppressWarnings(
    ComplexUpset::upset(
      mat,
      intersect = rev(methods),
      name = "Detection Method Overlap",
      min_size = min_size,
      sort_sets = FALSE,
      sort_intersections_by = "cardinality",
      height_ratio = 0.9,
      set_sizes = FALSE,
      base_annotations = list(
        "Detected Genes" =
          ComplexUpset::intersection_size(
            mapping = ggplot2::aes(fill = "bars_color"),
            counts = TRUE,
            text = list(size = 3, vjust = -0.35, hjust = -0.05, angle = 30)
          ) +
            ggplot2::scale_fill_manual(
              values = c("bars_color" = bar_col),
              guide = "none"
            )
      ),
      matrix = ComplexUpset::intersection_matrix(
        geom = ggplot2::geom_point(size = 3),
        segment = ggplot2::geom_segment(linewidth = 1.25),
        outline_color = list(
          active = bar_col,
          inactive = inactive_dot_col
        )
      ),
      stripes = ComplexUpset::upset_stripes(
        mapping = ggplot2::aes(color = layer),
        data = set_metadata,
        colors = stripe_cols
      ),
      themes = ComplexUpset::upset_modify_themes(
        list(
          "Detected Genes" = ggplot2::theme(
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_text(size = base_size),
            panel.grid.major = ggplot2::element_line(
              color = unname(grey_palette_colors[["light"]])
            ),
            panel.grid.minor = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(fill = "white", colour = NA),
            plot.background = ggplot2::element_rect(fill = "white", colour = NA)
          )
        )
      )
    ) +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(size = 18, face = "bold"),
        text = ggplot2::element_text(size = base_size)
      )
  )
}

plot_harmonized_set_sizes <- function(
  plot_data,
  timepoint,
  fdr = max(plot_data$meta$thresholds),
  base_size = 12
) {
  .check_pkg(c("dplyr", "ggplot2", "scales", "ggh4x"))

  sets <- dplyr::bind_rows(
    plot_data$threshold_sets$gene |>
      dplyr::filter(timepoint == !!timepoint, fdr == !!fdr),
    plot_data$threshold_sets$splicing |>
      dplyr::filter(timepoint == !!timepoint, fdr == !!fdr)
  )

  gene_genes <- sets |>
    dplyr::filter(layer == "Gene") |>
    dplyr::pull(ensgene) |>
    unique()

  splice_genes <- sets |>
    dplyr::filter(layer == "Splicing") |>
    dplyr::pull(ensgene) |>
    unique()

  overlap_genes <- intersect(gene_genes, splice_genes)

  method_df <- sets |>
    dplyr::group_by(method, layer) |>
    dplyr::summarise(n_genes = dplyr::n_distinct(ensgene), .groups = "drop") |>
    dplyr::mutate(
      panel = dplyr::case_when(
        layer == "Gene" ~ "Gene Level",
        layer == "Splicing" ~ "Transcript Level",
        TRUE ~ layer
      ),
      set = method
    ) |>
    dplyr::select(panel, set, n_genes)

  overlap_df <- data.frame(
    panel = "Overlap",
    set = .overlap_category_order(),
    n_genes = c(
      length(splice_genes),
      length(overlap_genes),
      length(gene_genes)
    ),
    stringsAsFactors = FALSE
  )

  df <- dplyr::bind_rows(method_df, overlap_df) |>
    dplyr::mutate(
      panel = factor(
        panel,
        levels = .overlap_category_order()
      ),
      set = factor(
        set,
        levels = unique(c(
          .splicing_methods(),
          .overlap_category_order(),
          .gene_methods()
        ))
      )
    )

  tp_cols <- .timepoint_shades(timepoint)

  fill_cols <- c(
    .tool_colors(shade = "base"),
    tp_cols[c("Transcript Level", "Overlap", "Gene Level")]
  )

  strip_cols <- tp_cols[.overlap_category_order()]

  fig <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = set, y = n_genes, fill = set)
  ) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::comma(n_genes)),
      vjust = -0.35,
      size = 3
    ) +
    ggh4x::facet_grid2(
      ~panel,
      scales = "free_x",
      space = "free_x",
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = unname(strip_cols[levels(df$panel)]),
          color = unname(grey_palette_colors[["black"]])
        ),
        text_x = ggh4x::elem_list_text(
          color = "white",
          face = "bold"
        )
      )
    ) +
    ggplot2::scale_x_discrete(drop = TRUE) +
    ggplot2::scale_fill_manual(values = fill_cols, drop = FALSE) +
    ggplot2::scale_y_continuous(
      labels = scales::comma,
      expand = ggplot2::expansion(mult = c(0, 0.18))
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Detected Genes"
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold")
    )

  fig
}

plot_threshold_overlap_bar <- function(
  plot_data,
  timepoint,
  base_size = 12
) {
  .check_pkg(c("dplyr", "ggplot2", "scales", "ggh4x"))

  df <- plot_data$threshold_overlap |>
    dplyr::filter(as.character(timepoint) == !!timepoint) |>
    dplyr::mutate(
      class = dplyr::case_when(
        class == "Gene only" ~ "Gene Level",
        class == "Splicing only" ~ "Transcript Level",
        class == "Overlap" ~ "Overlap",
        TRUE ~ as.character(class)
      ),
      class = factor(class, levels = .overlap_category_order()),
      panel = factor("FDR Threshold", levels = "FDR Threshold")
    )

  class_cols <- .timepoint_shades(timepoint)
  class_cols <- class_cols[.overlap_category_order()]

  tp_col <- .timepoint_main_color(timepoint)

  ggplot2::ggplot(df, ggplot2::aes(x = fdr_label, y = prop, fill = class)) +
    ggplot2::geom_col(width = 0.75) +
    ggh4x::facet_grid2(
      ~panel,
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = .timepoint_main_color(timepoint),
          color = unname(grey_palette_colors[["black"]])
        ),
        text_x = ggh4x::elem_list_text(
          color = "white",
          face = "bold"
        )
      )
    ) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::scale_fill_manual(
      values = class_cols,
      breaks = .overlap_category_order(),
      drop = FALSE,
      name = NULL
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Detected Genes (%)",
      fill = NULL
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none",
      legend.direction = "vertical",
      legend.key.size = grid::unit(0.45, "cm")
    )
}

plot_timepoint_row_title <- function(
  timepoint,
  base_size = 18,
  fill = "white"
) {
  .check_pkg("ggplot2")

  label <- .timepoint_labels()[[timepoint]]
  tp_col <- .timepoint_main_color(timepoint)

  ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 0.5,
      label = label,
      hjust = 0,
      fontface = "bold",
      color = tp_col,
      size = base_size / 3
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = fill, color = NA),
      plot.margin = ggplot2::margin(2, 2, 0, 2)
    )
}

plot_timepoint_shade_legend_one <- function(
  timepoint,
  base_size = 15,
  title_size = 18,
  key_size_cm = 0.75
) {
  .check_pkg(c("ggplot2", "cowplot"))

  cats <- .overlap_category_order()
  cols <- .timepoint_shades(timepoint)[cats]

  legend_df <- data.frame(
    category = factor(cats, levels = cats),
    x = 1,
    y = 1,
    stringsAsFactors = FALSE
  )

  legend_plot <- ggplot2::ggplot(
    legend_df,
    ggplot2::aes(x = x, y = y, fill = category)
  ) +
    ggplot2::geom_tile(alpha = 0) +
    ggplot2::scale_fill_manual(
      values = cols,
      name = .timepoint_labels()[[timepoint]],
      guide = ggplot2::guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        label.position = "right",
        nrow = 1,
        byrow = TRUE,
        override.aes = list(alpha = 1)
      )
    ) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_text(
        size = title_size,
        hjust = 0.5
      ),
      legend.text = ggplot2::element_text(size = base_size),
      legend.key.width = grid::unit(key_size_cm, "cm"),
      legend.key.height = grid::unit(key_size_cm, "cm"),
      legend.spacing.x = grid::unit(0.25, "cm"),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )

  cowplot::get_legend(legend_plot)
}

plot_timepoint_shade_legend <- function(
  base_size = 15,
  title_size = 18,
  key_size_cm = 0.75
) {
  .check_pkg(c("cowplot"))

  legends <- lapply(.timepoint_order(), function(tp) {
    plot_timepoint_shade_legend_one(
      timepoint = tp,
      base_size = base_size,
      title_size = title_size,
      key_size_cm = key_size_cm
    )
  })

  cowplot::plot_grid(
    plotlist = legends,
    nrow = 1,
    align = "h"
  )
}

plot_method_legend_one <- function(
  methods,
  colors,
  title,
  base_size = 15,
  title_size = 18,
  key_size_cm = 0.75
) {
  .check_pkg(c("ggplot2", "cowplot"))

  legend_df <- data.frame(
    method = factor(methods, levels = methods),
    x = 1,
    y = 1,
    stringsAsFactors = FALSE
  )

  legend_plot <- ggplot2::ggplot(
    legend_df,
    ggplot2::aes(x = x, y = y, fill = method)
  ) +
    ggplot2::geom_tile(alpha = 0) +
    ggplot2::scale_fill_manual(
      values = colors[methods],
      name = title,
      guide = ggplot2::guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        label.position = "right",
        nrow = 1,
        byrow = TRUE,
        override.aes = list(alpha = 1)
      )
    ) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_text(
        size = title_size,
        hjust = 0.5
      ),
      legend.text = ggplot2::element_text(size = base_size),
      legend.key.width = grid::unit(key_size_cm, "cm"),
      legend.key.height = grid::unit(key_size_cm, "cm"),
      legend.spacing.x = grid::unit(0.25, "cm"),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )

  cowplot::get_legend(legend_plot)
}

plot_method_legend <- function(
  base_size = 15,
  title_size = 18,
  key_size_cm = 0.75
) {
  .check_pkg("cowplot")

  cols <- .tool_colors(shade = "base")

  gene_legend <- plot_method_legend_one(
    methods = .gene_methods(),
    colors = cols,
    title = "Gene Level Methods",
    base_size = base_size,
    title_size = title_size,
    key_size_cm = key_size_cm
  )

  transcript_legend <- plot_method_legend_one(
    methods = .splicing_methods(),
    colors = cols,
    title = "Transcript Level Methods",
    base_size = base_size,
    title_size = title_size,
    key_size_cm = key_size_cm
  )

  cowplot::plot_grid(
    gene_legend,
    transcript_legend,
    nrow = 1,
    align = "h",
    rel_widths = c(1, 1)
  )
}

plot_filter_threshold_robustness_timepoint_panel <- function(
  plot_data,
  timepoint,
  upset_width = 4.25,
  set_size_width = 2.75,
  threshold_width = 1.5,
  min_size = 5,
  base_size = 12
) {
  .check_pkg(c("patchwork", "ggplot2"))

  tp_min_size <- .get_timepoint_min_size(min_size, timepoint)

  upset_plot <- plot_harmonized_method_upset(
    plot_data = plot_data,
    timepoint = timepoint,
    min_size = tp_min_size,
    base_size = base_size
  )

  set_size_plot <- plot_harmonized_set_sizes(
    plot_data = plot_data,
    timepoint = timepoint,
    base_size = base_size
  )

  threshold_plot <- plot_threshold_overlap_bar(
    plot_data = plot_data,
    timepoint = timepoint,
    base_size = base_size
  )

  panel <- patchwork::wrap_elements(full = upset_plot) +
    patchwork::wrap_elements(full = set_size_plot) +
    patchwork::wrap_elements(full = threshold_plot) +
    patchwork::plot_layout(
      widths = c(upset_width, set_size_width, threshold_width)
    )

  row_title <- plot_timepoint_row_title(timepoint)

  full_panel <- row_title / panel +
    patchwork::plot_layout(heights = c(0.25, 5))

  invisible(list(
    panel = full_panel,
    body = panel,
    upset = upset_plot,
    set_sizes = set_size_plot,
    threshold = threshold_plot
  ))
}

# --------------------------------------------------
# Convenience Wrapper
# --------------------------------------------------

prepare_splicing_robustness_inputs <- function(
  drimseq_results,
  timepoints = c("H1", "H3", "H24"),
  ref_level = "C1",
  rerun_dexseq = TRUE,
  prep_suppa_inputs = TRUE,
  load_suppa_results = FALSE,
  dexseq_force = FALSE,
  suppa_force = TRUE
) {
  message("[Sup Fig 4 Prep] Starting")

  # Harmonized universe
  message("[1/3] Building harmonized gene universe from DRIMSeq")

  universe_df <- .get_harmonized_universe_from_drimseq(drimseq_results)

  gene_universe_by_tp <- split(
    universe_df$ensgene,
    universe_df$timepoint
  )

  universe_sizes <- vapply(gene_universe_by_tp, length, integer(1))

  message("[Universe] Genes per timepoint:")
  for (tp in names(universe_sizes)) {
    message("  - ", tp, ": ", universe_sizes[[tp]])
  }

  # DEXSeq harmonized
  message("\n[2/3] DEXSeq harmonized results...")

  dexseq_harmonized_results <- NULL
  dexseq_rds <- file.path(
    .get_results_dir(),
    "analysis",
    "dexseq",
    "dexseq_results.C1.summarizeOverlaps.multiOverlap.all.harmonized.rds"
  )

  if (isTRUE(rerun_dexseq)) {
    message("  -> Running DEXSeq on harmonized universe...")
    message("     Output: ", dexseq_rds)

    dexseq_harmonized_results <- run_dexseq(
      timepoints = timepoints,
      gene_universe_by_tp = gene_universe_by_tp,
      force = dexseq_force,
      out_rds = dexseq_rds,
      out_xlsx = sub("\\.rds$", ".xlsx", dexseq_rds)
    )
  } else if (file.exists(dexseq_rds)) {
    message("  -> Loading existing harmonized DEXSeq results...")
    message("     ", dexseq_rds)

    dexseq_harmonized_results <- readRDS(dexseq_rds)
  } else {
    stop("Harmonized DEXSeq RDS not found: ", dexseq_rds)
  }

  # SUPPA2 prep
  message("\n[3/3] SUPPA2 harmonized input preparation...")

  suppa_prep <- NULL

  if (isTRUE(prep_suppa_inputs)) {
    message("  -> Generating filtered isoTPM inputs...")

    suppa_prep <- .prepare_suppa_harmonized_isotpm(
      drimseq_results = drimseq_results
    )

    message("\n  [SUPPA2] Filtered isoTPM summary:")
    print(suppa_prep$isotpm_paths[, c("timepoint", "n_ref_tx", "n_test_tx")])

    message(
      "\n[SUPPA2] Harmonized isoTPM inputs written to:\n  ",
      suppa_prep$outdir,
      "\n\nRun SUPPA2 from the terminal with:\n\n",
      '"$PROJECT_ROOT/scripts/run" bash ',
      '"$PROJECT_ROOT/scripts/analysis/run_suppa.sh" ',
      ref_level, " harmonized ",
      '"', suppa_prep$outdir, '"',
      "\n"
    )
  } else {
    message("  -> Skipping isoTPM generation (prep_suppa_inputs = FALSE)")
  }

  # SUPPA2 load
  suppa_harmonized_results <- NULL

  if (isTRUE(load_suppa_results)) {
    suppa_dir <- file.path(
      .get_results_dir(),
      "analysis",
      "suppa",
      "harmonized",
      ref_level
    )

    message("\n[SUPPA2] Loading harmonized results...")
    message("  -> Directory: ", suppa_dir)

    suppa_harmonized_results <- load_suppa(
      ref_level = ref_level,
      suppa_dir = suppa_dir,
      force = suppa_force
    )

    message("  -> Load complete")
  } else {
    message("\n[SUPPA2] Skipping result load (load_suppa_results = FALSE)")
  }

  message("\n[Sup Fig 4 Prep] Complete\n")

  list(
    universe_df = universe_df,
    gene_universe_by_tp = gene_universe_by_tp,
    dexseq_harmonized_results = dexseq_harmonized_results,
    suppa_prep = suppa_prep,
    suppa_harmonized_results = suppa_harmonized_results
  )
}

plot_filter_threshold_robustness_all <- function(
  dexseq_results = NULL,
  drimseq_results = NULL,
  suppa_results = NULL,
  plot_data = NULL,
  outdir = NULL,
  plot_data_rds = NULL,
  use_cached_plot_data = TRUE,
  save_plot_data = TRUE,
  regenerate_plots = TRUE,
  min_size = c(H1 = 1, H3 = 5, H24 = 5),
  thresholds = c(0.01, 0.05, 0.10),
  save_pdf = TRUE,
  filename = "sup_fig_4_filter_threshold_robustness_analysis.pdf",
  ...
) {
  .check_pkg(c("patchwork", "ggplot2", "cowplot"))

  outdir <- .set_outdir(outdir, subdir = "sup_fig_4")

  if (is.null(plot_data_rds)) {
    plot_data_rds <- file.path(outdir, "sup_fig_4_filter_threshold_robustness_analysis.rds")
  }

  if (is.null(plot_data)) {
    if (isTRUE(use_cached_plot_data) && file.exists(plot_data_rds)) {
      message("Loading cached plot data: ", plot_data_rds)
      plot_data <- readRDS(plot_data_rds)
    } else {
      plot_data <- make_filter_threshold_robustness_data(
        dexseq_results = dexseq_results,
        drimseq_results = drimseq_results,
        suppa_results = suppa_results,
        thresholds = thresholds,
        ...
      )

      if (isTRUE(save_plot_data)) {
        saveRDS(plot_data, plot_data_rds)
        message("Saved plot data: ", plot_data_rds)
      }
    }
  }

  figs <- NULL
  if (isTRUE(regenerate_plots)) {
    timepoint_panels <- lapply(.timepoint_order(), function(tp) {
      plot_filter_threshold_robustness_timepoint_panel(
        plot_data = plot_data,
        timepoint = tp,
        min_size = min_size
      )
    })
    names(timepoint_panels) <- .timepoint_order()

    row_heights <- rep(1, length(timepoint_panels))

    combined_body <- patchwork::wrap_plots(
      lapply(timepoint_panels, `[[`, "panel"),
      ncol = 1,
      heights = row_heights
    )

    timepoint_legend <- plot_timepoint_shade_legend(
      base_size = 15,
      title_size = 18,
      key_size_cm = 0.75
    )

    method_legend <- plot_method_legend(
      base_size = 15,
      title_size = 18,
      key_size_cm = 0.75
    )

    legend_panel <- patchwork::wrap_plots(
      patchwork::wrap_elements(full = method_legend),
      patchwork::wrap_elements(full = timepoint_legend),
      ncol = 1,
      heights = c(1, 1)
    )

    body_height <- 5.5 * length(timepoint_panels)
    legend_height <- 2

    combined <- patchwork::wrap_plots(
      patchwork::wrap_elements(full = combined_body),
      patchwork::wrap_elements(full = legend_panel),
      ncol = 1,
      heights = c(body_height, legend_height)
    )

    ggplot2::ggsave(
      filename = file.path(outdir, filename),
      plot = combined,
      width = 15.5,
      height = body_height + legend_height
    )

    figs <- list(
      timepoint_panels = timepoint_panels,
      combined = combined
    )
  }

  invisible(list(
    plot_data = plot_data,
    figs = figs,
    paths = list(
      outdir = outdir,
      plot_data_rds = plot_data_rds,
      pdf = file.path(outdir, filename)
    )
  ))
}
