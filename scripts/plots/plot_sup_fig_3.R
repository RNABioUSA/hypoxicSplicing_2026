# ============================================================
# plot_sup_fig_3.R
# ============================================================

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

.drimseq_timepoint_order <- function() {
  c("H1", "H3", "H24")
}

.drimseq_timepoint_labels <- function() {
  stats::setNames(
    names(timepoint_base_colors)[grepl("Hypoxia", names(timepoint_base_colors))],
    .drimseq_timepoint_order()
  )
}

.drimseq_timepoint_colors <- function(shade = "light") {
  timepoints <- .drimseq_timepoint_order()

  tp_shades <- .get_timepoint_shaded_colors(
    timepoints = timepoints,
    labels = c("light", "base", "dark")
  )

  if (!shade %in% names(tp_shades[[1]])) {
    stop(
      "Invalid DRIMSeq timepoint shade: ", shade,
      ". Valid options are: ", paste(names(tp_shades[[1]]), collapse = ", ")
    )
  }

  stats::setNames(
    vapply(tp_shades, function(x) unname(x[[shade]]), character(1)),
    timepoints
  )
}

.drimseq_timepoint_strip <- function(shade = "light") {
  .check_pkg("ggh4x")

  ggh4x::strip_themed(
    background_x = ggh4x::elem_list_rect(
      fill = unname(.drimseq_timepoint_colors(shade = shade)),
      color = unname(grey_palette_colors["black"])
    ),
    text_x = ggh4x::elem_list_text(
      color = unname(grey_palette_colors["black"]),
      face = "bold"
    )
  )
}

.drimseq_tx_colors <- function(n) {
  category_ids <- paste("Category", 3:5)

  if (n > length(category_ids)) {
    warning(
      "Requested ", n, " transcript colors but only ", length(category_ids),
      " selected category colors are available. Colors will be recycled."
    )
  }

  category_shades <- .generate_shaded_palette(
    base_colors = category_base_colors[category_ids],
    labels = c("light", "base", "dark")
  )

  base_cols <- vapply(
    category_shades,
    function(x) unname(x["light"]),
    character(1)
  )

  cols <- rep(base_cols, length.out = n)
  stats::setNames(cols, paste0("Tx ", seq_len(n)))
}

.drimseq_gene_sets_by_tp <- function(drimseq_results) {
  timepoints <- .drimseq_timepoint_order()

  genes_by_tp <- lapply(timepoints, function(tp) {
    if (!is.null(drimseq_results$results$gene_sig_by_tp[[tp]])) {
      unique(stats::na.omit(drimseq_results$results$gene_sig_by_tp[[tp]]$ensgene))
    } else if (!is.null(drimseq_results$by_tp[[tp]]$results$gene_sig)) {
      unique(stats::na.omit(drimseq_results$by_tp[[tp]]$results$gene_sig$ensgene))
    } else if (!is.null(drimseq_results$by_tp[[tp]]$gene_sig)) {
      unique(stats::na.omit(drimseq_results$by_tp[[tp]]$gene_sig$ensgene))
    } else {
      character(0)
    }
  })
  names(genes_by_tp) <- timepoints
  genes_by_tp
}

.extract_drimseq_proportions_by_tp <- function(drimseq_results, ref_level = NULL) {
  .check_pkg(c("dplyr", "tidyr", "purrr", "DRIMSeq"))

  if (!("by_tp" %in% names(drimseq_results))) {
    stop("Expected drimseq_results$by_tp.")
  }

  ref_level <- ref_level %||% drimseq_results$meta$ref_level
  if (is.null(ref_level)) stop("Could not infer DRIMSeq reference level.")

  purrr::imap_dfr(drimseq_results$by_tp, function(tp_obj, tp) {
    if (is.null(tp_obj$dm)) stop("Expected drimseq_results$by_tp[[", tp, "]]$dm.")
    if (is.null(tp_obj$results$tx_full)) {
      stop("Expected drimseq_results$by_tp[[", tp, "]]$results$tx_full.")
    }

    d <- tp_obj$dm
    tx_full <- as.data.frame(tp_obj$results$tx_full, stringsAsFactors = FALSE)

    prop <- DRIMSeq::proportions(d)
    prop_df <- as.data.frame(prop, check.names = FALSE, stringsAsFactors = FALSE)

    samp <- as.data.frame(DRIMSeq::samples(d), stringsAsFactors = FALSE)
    sample_col <- c("sample_id", "sample", "names")[c("sample_id", "sample", "names") %in% names(samp)][1]

    if (is.na(sample_col)) stop("Could not identify sample column in DRIMSeq::samples(d).")
    if (!("condition" %in% names(samp))) stop("Could not identify 'condition' column in DRIMSeq::samples(d).")

    sample_ids <- as.character(samp[[sample_col]])

    if (!("feature_id" %in% names(prop_df))) {
      prop_df$feature_id <- rownames(prop_df)
    }

    prop_sample_cols <- intersect(sample_ids, names(prop_df))
    if (length(prop_sample_cols) == 0) {
      stop("No sample columns from DRIMSeq::samples(d) were found in DRIMSeq::proportions(d).")
    }

    long <- prop_df |>
      tidyr::pivot_longer(
        cols = tidyselect::all_of(prop_sample_cols),
        names_to = "sample_id",
        values_to = "proportion"
      ) |>
      dplyr::left_join(
        samp |>
          dplyr::transmute(sample_id = as.character(.data[[sample_col]]), condition = as.character(condition)),
        by = "sample_id"
      ) |>
      dplyr::mutate(timepoint = tp)

    tx_meta <- tx_full |>
      dplyr::distinct(feature_id, ensgene, symbol, padj, significant, .keep_all = TRUE) |>
      dplyr::select(feature_id, ensgene, symbol, padj, significant)

    long |>
      dplyr::left_join(tx_meta, by = "feature_id") |>
      dplyr::mutate(
        symbol = dplyr::if_else(is.na(symbol) | symbol == "", ensgene, symbol),
        comparison_group = dplyr::case_when(
          condition == ref_level ~ ref_level,
          TRUE ~ tp
        ),
        timepoint = factor(timepoint, levels = .drimseq_timepoint_order())
      )
  })
}

.summarize_drimseq_mean_proportions <- function(proportion_df) {
  .check_pkg("dplyr")

  proportion_df |>
    dplyr::group_by(timepoint, comparison_group, ensgene, symbol, feature_id, padj, significant) |>
    dplyr::summarize(
      mean_prop = mean(proportion, na.rm = TRUE),
      .groups = "drop"
    )
}

.summarize_drimseq_delta_proportions <- function(mean_prop_df, ref_level) {
  .check_pkg("dplyr")

  ref_df <- mean_prop_df |>
    dplyr::filter(comparison_group == ref_level) |>
    dplyr::select(timepoint, ensgene, feature_id, ref_prop = mean_prop)

  test_df <- mean_prop_df |>
    dplyr::filter(comparison_group != ref_level) |>
    dplyr::select(
      timepoint, comparison_group, ensgene, symbol,
      feature_id, padj, significant,
      test_prop = mean_prop
    )

  test_df |>
    dplyr::left_join(ref_df, by = c("timepoint", "ensgene", "feature_id")) |>
    dplyr::mutate(
      dprop = test_prop - ref_prop,
      abs_dprop = abs(dprop)
    )
}

.make_drimseq_multiplicity_magnitude_df <- function(delta_prop_df, sig_only = TRUE) {
  .check_pkg("dplyr")

  dprop_df <- delta_prop_df

  if (isTRUE(sig_only)) {
    dprop_df <- dprop_df |>
      dplyr::filter(significant %in% TRUE)
  }

  dprop_df |>
    dplyr::filter(
      !is.na(ensgene),
      !is.na(feature_id),
      is.finite(abs_dprop)
    ) |>
    dplyr::group_by(timepoint, ensgene, symbol) |>
    dplyr::summarize(
      n_sig_tx = dplyr::n_distinct(feature_id),
      max_abs_dprop = max(abs_dprop, na.rm = TRUE),
      mean_abs_dprop = mean(abs_dprop, na.rm = TRUE),
      min_padj = suppressWarnings(min(padj, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      min_padj = dplyr::if_else(is.infinite(min_padj), NA_real_, min_padj),
      tool = "DRIMSeq"
    )
}

.select_drimseq_top_genes_by_timepoint <- function(
  delta_prop_df,
  n_per_timepoint = 5
) {
  .check_pkg("dplyr")

  delta_prop_df |>
    dplyr::filter(significant %in% TRUE, is.finite(abs_dprop)) |>
    dplyr::group_by(timepoint, ensgene, symbol) |>
    dplyr::summarize(
      n_sig_tx = dplyr::n_distinct(feature_id),
      max_abs_dprop = max(abs_dprop, na.rm = TRUE),
      mean_abs_dprop = mean(abs_dprop, na.rm = TRUE),
      min_padj = suppressWarnings(min(padj, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      min_padj = dplyr::if_else(is.infinite(min_padj), NA_real_, min_padj)
    ) |>
    dplyr::group_by(timepoint) |>
    dplyr::arrange(
      dplyr::desc(max_abs_dprop),
      dplyr::desc(n_sig_tx),
      min_padj,
      .by_group = TRUE
    ) |>
    dplyr::slice_head(n = n_per_timepoint) |>
    dplyr::ungroup()
}

.make_drimseq_top_transcripts_df <- function(
  delta_prop_df,
  top_genes_df,
  top_tx_per_gene = 3
) {
  .check_pkg("dplyr")

  top_genes_df <- top_genes_df |>
    dplyr::rename(selected_timepoint = timepoint)

  delta_prop_df |>
    dplyr::inner_join(
      top_genes_df |>
        dplyr::select(selected_timepoint, ensgene, symbol),
      by = c("ensgene", "symbol"),
      relationship = "many-to-many"
    ) |>
    dplyr::filter(timepoint == selected_timepoint) |>
    dplyr::group_by(selected_timepoint, ensgene, symbol, feature_id) |>
    dplyr::summarize(
      max_abs_dprop = max(abs_dprop, na.rm = TRUE),
      mean_abs_dprop = mean(abs_dprop, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::group_by(selected_timepoint, ensgene, symbol) |>
    dplyr::arrange(
      dplyr::desc(max_abs_dprop),
      dplyr::desc(mean_abs_dprop),
      .by_group = TRUE
    ) |>
    dplyr::slice_head(n = top_tx_per_gene) |>
    dplyr::ungroup()
}

.make_drimseq_trajectory_df <- function(
  common_mean_prop_df,
  top_tx_df,
  condition_order = c("C1", "H1", "H3", "H24")
) {
  .check_pkg("dplyr")

  selected_tx <- top_tx_df |>
    dplyr::distinct(
      selected_timepoint,
      ensgene,
      feature_id,
      symbol
    ) |>
    dplyr::rename(selected_symbol = symbol)

  out <- common_mean_prop_df |>
    dplyr::select(-dplyr::any_of("symbol")) |>
    dplyr::inner_join(
      selected_tx,
      by = c("ensgene", "feature_id"),
      relationship = "many-to-many"
    ) |>
    dplyr::mutate(
      symbol = selected_symbol,
      condition = factor(
        as.character(condition),
        levels = condition_order
      ),
      selected_timepoint = factor(
        as.character(selected_timepoint),
        levels = .drimseq_timepoint_order(),
        labels = unname(.drimseq_timepoint_labels())
      ),
      comparison_group = factor(
        as.character(condition),
        levels = condition_order
      )
    ) |>
    dplyr::group_by(
      selected_timepoint,
      ensgene,
      symbol
    ) |>
    dplyr::mutate(
      tx_rank = match(
        feature_id,
        sort(unique(feature_id))
      ),
      tx_label = paste0("Tx ", tx_rank)
    ) |>
    dplyr::ungroup()

  coverage_check <- out |>
    dplyr::distinct(
      selected_timepoint,
      ensgene,
      feature_id,
      condition
    ) |>
    dplyr::count(
      selected_timepoint,
      ensgene,
      feature_id,
      name = "n_conditions"
    ) |>
    dplyr::filter(n_conditions != length(condition_order))

  if (nrow(coverage_check) > 0) {
    stop(
      nrow(coverage_check),
      " selected transcript trajectories do not contain all ",
      length(condition_order),
      " conditions."
    )
  }

  out
}

.make_drimseq_trajectory_label_df <- function(trajectory_df) {
  .check_pkg("dplyr")

  trajectory_df |>
    dplyr::group_by(selected_timepoint, ensgene, symbol, feature_id) |>
    dplyr::arrange(comparison_group, .by_group = TRUE) |>
    dplyr::slice_tail(n = 1) |>
    dplyr::ungroup()
}

.reconstruct_drimseq_common_counts <- function(
  drimseq_results,
  condition_order = c("C1", "H1", "H3", "H24")
) {
  .check_pkg(c("dplyr", "tidyr", "purrr", "DRIMSeq"))

  if (is.null(drimseq_results$by_tp)) {
    stop("Expected drimseq_results$by_tp.")
  }

  extracted <- purrr::imap(
    drimseq_results$by_tp,
    function(tp_obj, tp) {
      if (is.null(tp_obj$counts)) {
        stop("Missing stored counts for DRIMSeq comparison: ", tp)
      }

      if (is.null(tp_obj$dm)) {
        stop("Missing stored dm object for DRIMSeq comparison: ", tp)
      }

      counts_df <- as.data.frame(
        tp_obj$counts,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )

      samples_df <- as.data.frame(
        DRIMSeq::samples(tp_obj$dm),
        stringsAsFactors = FALSE
      )

      sample_col <- c("sample_id", "sample", "names")[
        c("sample_id", "sample", "names") %in% names(samples_df)
      ][1]

      if (is.na(sample_col)) {
        stop("Could not identify the DRIMSeq sample column for ", tp)
      }

      if (!("condition" %in% names(samples_df))) {
        stop("Could not identify the condition column for ", tp)
      }

      sample_ids <- as.character(samples_df[[sample_col]])
      count_sample_cols <- intersect(sample_ids, names(counts_df))

      if (length(count_sample_cols) == 0) {
        stop("No stored count columns matched samples for ", tp)
      }

      count_long <- counts_df |>
        dplyr::select(
          gene_id,
          feature_id,
          dplyr::all_of(count_sample_cols)
        ) |>
        tidyr::pivot_longer(
          cols = dplyr::all_of(count_sample_cols),
          names_to = "sample_id",
          values_to = "count"
        ) |>
        dplyr::mutate(
          gene_id = as.character(gene_id),
          feature_id = as.character(feature_id),
          sample_id = as.character(sample_id),
          source_comparison = tp
        )

      sample_meta <- samples_df |>
        dplyr::transmute(
          sample_id = as.character(.data[[sample_col]]),
          condition = as.character(condition)
        )

      list(
        counts = count_long,
        samples = sample_meta
      )
    }
  )

  sample_df <- dplyr::bind_rows(
    purrr::map(extracted, "samples")
  ) |>
    dplyr::distinct(sample_id, condition) |>
    dplyr::mutate(
      condition = factor(condition, levels = condition_order)
    )

  # Normoxic counts occur in all three pairwise tables.
  # They should be identical, so collapse duplicated entries.
  duplicated_count_check <- dplyr::bind_rows(
    purrr::map(extracted, "counts")
  ) |>
    dplyr::group_by(gene_id, feature_id, sample_id) |>
    dplyr::summarize(
      n_values = dplyr::n_distinct(count),
      .groups = "drop"
    ) |>
    dplyr::filter(n_values > 1)

  if (nrow(duplicated_count_check) > 0) {
    stop(
      nrow(duplicated_count_check),
      " transcript/sample entries differed between cached pairwise count tables."
    )
  }

  observed_counts <- dplyr::bind_rows(
    purrr::map(extracted, "counts")
  ) |>
    dplyr::group_by(gene_id, feature_id, sample_id) |>
    dplyr::summarize(
      count = dplyr::first(count),
      .groups = "drop"
    )

  transcript_map <- observed_counts |>
    dplyr::distinct(gene_id, feature_id)

  # Complete every retained transcript across every sample.
  # Missing transcript/sample combinations represent zero counts.
  common_counts_long <- tidyr::crossing(
    transcript_map,
    sample_id = sample_df$sample_id
  ) |>
    dplyr::left_join(
      observed_counts,
      by = c("gene_id", "feature_id", "sample_id")
    ) |>
    dplyr::mutate(
      count = dplyr::coalesce(count, 0)
    ) |>
    dplyr::left_join(sample_df, by = "sample_id")

  list(
    counts_long = common_counts_long,
    sample_df = sample_df,
    duplicate_check = duplicated_count_check
  )
}

.make_common_drimseq_mean_proportions <- function(
  drimseq_results,
  annot_df = NULL,
  condition_order = c("C1", "H1", "H3", "H24")
) {
  .check_pkg(c("dplyr", "tidyr"))

  common <- .reconstruct_drimseq_common_counts(
    drimseq_results = drimseq_results,
    condition_order = condition_order
  )

  proportions <- common$counts_long |>
    dplyr::group_by(sample_id, condition, gene_id) |>
    dplyr::mutate(
      gene_count = sum(count, na.rm = TRUE),
      proportion = dplyr::if_else(
        gene_count > 0,
        count / gene_count,
        0
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      ensgene = .strip_ens_version(gene_id)
    )

  if (!is.null(annot_df)) {
    gene_annot <- annot_df |>
      dplyr::distinct(ensgene, .keep_all = TRUE) |>
      dplyr::select(ensgene, symbol)

    proportions <- proportions |>
      dplyr::left_join(gene_annot, by = "ensgene")
  }

  if (!("symbol" %in% names(proportions))) {
    proportions$symbol <- proportions$ensgene
  }

  proportions |>
    dplyr::mutate(
      symbol = dplyr::if_else(
        is.na(symbol) | symbol == "",
        ensgene,
        symbol
      )
    ) |>
    dplyr::group_by(
      condition,
      ensgene,
      symbol,
      feature_id
    ) |>
    dplyr::summarize(
      mean_prop = mean(proportion, na.rm = TRUE),
      n_samples = dplyr::n_distinct(sample_id),
      .groups = "drop"
    )
}

# --------------------------------------------------
# Main Analysis
# --------------------------------------------------

run_drimseq_supplement_analysis <- function(
  drimseq_results,
  outdir = NULL,
  sig_only = TRUE,
  ref_level = NULL,
  n_genes_per_timepoint = 5,
  top_tx_per_gene = 3,
  save_tables = FALSE,
  force_recompute = FALSE
) {
  .check_pkg(c("dplyr", "tidyr", "purrr", "DRIMSeq"))

  outdir <- .set_outdir(outdir, subdir = "sup_fig_3")
  outfile <- file.path(outdir, "sup_fig_3_drimseq_analysis.rds")

  if (!isTRUE(force_recompute) && file.exists(outfile)) {
    message("Loading cached DRIMSeq supplement analysis from: ", outfile)
    return(readRDS(outfile))
  }

  ref_level <- ref_level %||% drimseq_results$meta$ref_level
  if (is.null(ref_level)) stop("Could not infer DRIMSeq reference level.")

  genes_by_tp <- .drimseq_gene_sets_by_tp(drimseq_results)

  proportion_df <- .extract_drimseq_proportions_by_tp(
    drimseq_results = drimseq_results,
    ref_level = ref_level
  )

  mean_prop_df <- .summarize_drimseq_mean_proportions(proportion_df)

  delta_prop_df <- .summarize_drimseq_delta_proportions(
    mean_prop_df = mean_prop_df,
    ref_level = ref_level
  )

  multiplicity_df <- .make_drimseq_multiplicity_magnitude_df(
    delta_prop_df = delta_prop_df,
    sig_only = sig_only
  )

  top_genes_df <- .select_drimseq_top_genes_by_timepoint(
    delta_prop_df = delta_prop_df,
    n_per_timepoint = n_genes_per_timepoint
  )

  top_tx_df <- .make_drimseq_top_transcripts_df(
    delta_prop_df = delta_prop_df,
    top_genes_df = top_genes_df,
    top_tx_per_gene = top_tx_per_gene
  )

  common_mean_prop_df <- .make_common_drimseq_mean_proportions(
    drimseq_results = drimseq_results,
    annot_df = .get_annot(),
    condition_order = c("C1", "H1", "H3", "H24")
  )

  trajectory_df <- .make_drimseq_trajectory_df(
    common_mean_prop_df = common_mean_prop_df,
    top_tx_df = top_tx_df,
    condition_order = c("C1", "H1", "H3", "H24")
  )

  trajectory_label_df <- .make_drimseq_trajectory_label_df(trajectory_df)

  out <- list(
    meta = list(
      outdir = outdir,
      sig_only = sig_only,
      ref_level = ref_level,
      n_genes_per_timepoint = n_genes_per_timepoint,
      top_tx_per_gene = top_tx_per_gene,
      timepoint_order = .drimseq_timepoint_order()
    ),
    overlap = list(
      genes_by_tp = genes_by_tp
    ),
    proportions = list(
      raw_proportion_df = proportion_df,
      mean_prop_df = mean_prop_df,
      delta_prop_df = delta_prop_df
    ),
    multiplicity = list(
      gene_df = multiplicity_df
    ),
    trajectories = list(
      top_genes_df = top_genes_df,
      top_tx_df = top_tx_df,
      plot_df = trajectory_df,
      label_df = trajectory_label_df
    )
  )

  if (isTRUE(save_tables)) {
    utils::write.table(
      multiplicity_df,
      file = file.path(outdir, "sup_fig_3_B_drimseq_multiplicity_magnitude.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    utils::write.table(
      top_genes_df,
      file = file.path(outdir, "sup_fig_3_C_drimseq_top_genes.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    utils::write.table(
      trajectory_df,
      file = file.path(outdir, "sup_fig_3_C_drimseq_transcript_trajectories.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  saveRDS(out, outfile)
  out
}

# --------------------------------------------------
# Figure S2A - DRIMSeq Gene Overlap Across Timepoints
# --------------------------------------------------

plot_drimseq_overlap_across_timepoints <- function(
  res,
  outdir = NULL,
  title = "DRIMSeq: Gene Overlap Across Hypoxia Timepoints"
) {
  .check_pkg("grid")

  outdir <- .set_outdir(outdir, subdir = "sup_fig_3")

  fig_a <- plot_feature_overlap_across_timepoints(
    features_by_tp = res$overlap$genes_by_tp,
    title = title,
    shade_label = "TX"
  )

  grid::grid.newpage()
  grid::grid.draw(fig_a$plot)

  grDevices::pdf(
    file.path(outdir, "sup_fig_3_A_drimseq_gene_overlap.pdf"),
    width = 6.75,
    height = 6.75
  )
  grid::grid.newpage()
  grid::grid.draw(fig_a$plot)
  invisible(grDevices::dev.off())

  invisible(fig_a)
}

# --------------------------------------------------
# Figure S2B - DRIMSeq Multiplicity and Magnitude
# --------------------------------------------------

plot_drimseq_multiplicity_vs_magnitude <- function(
  res,
  outdir = NULL,
  title = "DRIMSeq: Transcript Multiplicity and Magnitude",
  base_size = 15,
  max_features_to_show = 10,
  max_effect_to_show = 1,
  label_top_n_by_features = 5,
  label_min_features = 5,
  point_alpha = 0.75,
  point_size = 3
) {
  .check_pkg(c("dplyr", "ggplot2", "ggrepel", "ggh4x", "scales"))

  outdir <- .set_outdir(outdir, subdir = "sup_fig_3")

  df <- res$multiplicity$gene_df

  plot_df <- df |>
    dplyr::filter(
      n_sig_tx <= max_features_to_show,
      max_abs_dprop <= max_effect_to_show
    )

  excluded_df <- df |>
    dplyr::filter(
      n_sig_tx > max_features_to_show |
        max_abs_dprop > max_effect_to_show
    )

  global_median_features <- stats::median(plot_df$n_sig_tx, na.rm = TRUE)
  global_median_effect <- stats::median(plot_df$max_abs_dprop, na.rm = TRUE)

  global_subtitle <- bquote(
    "Median Transcripts per Gene = " * .(round(global_median_features, 1)) *
      "; Median Max |" * Delta * " Proportion| = " *
      .(round(global_median_effect, 2))
  )

  stat_df <- plot_df |>
    dplyr::group_by(timepoint) |>
    dplyr::summarize(
      median_features = stats::median(n_sig_tx, na.rm = TRUE),
      mean_features = mean(n_sig_tx, na.rm = TRUE),
      median_effect = stats::median(max_abs_dprop, na.rm = TRUE),
      mean_effect = mean(max_abs_dprop, na.rm = TRUE),
      n_genes = dplyr::n(),
      .groups = "drop"
    )

  stat_label_df <- stat_df |>
    dplyr::mutate(
      stat_label = paste0(
        "atop(",
        "'Mean Transcripts = ", round(mean_features, 2), "',",
        "'Mean |' * Delta * ' Proportion| = ", round(mean_effect, 3), "'",
        ")"
      )
    )

  label_df <- plot_df |>
    dplyr::group_by(timepoint) |>
    dplyr::filter(n_sig_tx >= label_min_features) |>
    dplyr::arrange(dplyr::desc(n_sig_tx), dplyr::desc(max_abs_dprop), .by_group = TRUE) |>
    dplyr::slice_head(n = label_top_n_by_features) |>
    dplyr::ungroup() |>
    dplyr::mutate(label_symbol = symbol)

  repel_df <- plot_df |>
    dplyr::left_join(
      label_df |>
        dplyr::select(timepoint, ensgene, label_symbol),
      by = c("timepoint", "ensgene")
    )

  fig_b <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = n_sig_tx, y = max_abs_dprop)
  ) +
    ggplot2::geom_vline(
      data = stat_df,
      ggplot2::aes(xintercept = median_features),
      linetype = "dashed",
      linewidth = 0.5,
      color = unname(grey_palette_colors["dark"])
    ) +
    ggplot2::geom_hline(
      data = stat_df,
      ggplot2::aes(yintercept = median_effect),
      linetype = "dashed",
      linewidth = 0.5,
      color = unname(grey_palette_colors["dark"])
    ) +
    ggplot2::geom_label(
      data = stat_label_df,
      ggplot2::aes(x = Inf, y = -Inf, label = stat_label),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = -0.35,
      size = 3,
      linewidth = 0.3,
      label.padding = grid::unit(0.18, "lines"),
      fill = "white",
      alpha = 0.9,
      parse = TRUE
    ) +
    ggplot2::geom_point(
      color = splicing_tool_colors[["DRIMSeq"]],
      alpha = point_alpha,
      size = point_size
    ) +
    ggrepel::geom_text_repel(
      data = repel_df,
      ggplot2::aes(label = label_symbol),
      size = 3,
      max.overlaps = Inf,
      max.iter = 50000,
      max.time = 5,
      force = 3,
      force_pull = 0.15,
      box.padding = 0.75,
      point.padding = 0.75,
      min.segment.length = 0.1,
      segment.alpha = 0.5,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    ggh4x::facet_wrap2(
      ~timepoint,
      nrow = 1,
      labeller = ggplot2::labeller(timepoint = .drimseq_timepoint_labels()),
      strip = .drimseq_timepoint_strip(shade = "light")
    ) +
    ggplot2::scale_x_continuous(
      limits = c(1, max_features_to_show),
      breaks = 1:max_features_to_show
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, max_effect_to_show),
      breaks = scales::pretty_breaks(n = 6)
    ) +
    ggplot2::labs(
      x = "Significant DRIMSeq Transcripts per Gene",
      y = expression("Max |" * Delta * " Transcript Proportion| per Gene"),
      title = title,
      subtitle = global_subtitle
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(
        linewidth = 0.5,
        color = unname(grey_palette_colors["light"])
      ),
      panel.grid.major.y = ggplot2::element_line(
        linewidth = 0.5,
        color = unname(grey_palette_colors["light"])
      ),
      axis.text.x = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 12)
    )

  attr(fig_b, "data_all") <- df
  attr(fig_b, "data_plotted") <- plot_df
  attr(fig_b, "excluded_genes") <- excluded_df
  attr(fig_b, "stats") <- stat_df

  ggplot2::ggsave(
    filename = file.path(outdir, "sup_fig_3_B_drimseq_multiplicity_magnitude.pdf"),
    plot = fig_b,
    width = 12.5,
    height = 9
  )

  fig_b
}

# --------------------------------------------------
# Figure S2C - DRIMSeq Transcript Usage Trajectories
# --------------------------------------------------

plot_drimseq_transcript_trajectories <- function(
  res,
  outdir = NULL,
  title = "DRIMSeq: Transcript Usage Trajectories",
  base_size = 15,
  rel_heights = c(1, 1, 1)
) {
  .check_pkg(c(
    "dplyr",
    "ggplot2",
    "ggrepel",
    "scales",
    "patchwork",
    "ggh4x"
  ))

  outdir <- .set_outdir(outdir, subdir = "sup_fig_3")

  plot_df <- res$trajectories$plot_df
  label_df <- res$trajectories$label_df

  n_genes_per_timepoint <- res$meta$n_genes_per_timepoint
  top_tx_per_gene <- res$meta$top_tx_per_gene

  tx_colors <- .drimseq_tx_colors(top_tx_per_gene)

  make_tp_plot <- function(tp_raw, show_title = FALSE) {
    tp_label <- .drimseq_timepoint_labels()[[tp_raw]]

    plot_tp <- plot_df |>
      dplyr::filter(as.character(selected_timepoint) == tp_label) |>
      dplyr::mutate(symbol = droplevels(factor(symbol)))

    label_tp <- label_df |>
      dplyr::filter(as.character(selected_timepoint) == tp_label) |>
      dplyr::mutate(symbol = droplevels(factor(symbol)))

    tp_color <- unname(.drimseq_timepoint_colors(shade = "light")[tp_raw])

    ggplot2::ggplot(
      plot_tp,
      ggplot2::aes(
        x = comparison_group,
        y = mean_prop,
        color = tx_label,
        group = feature_id
      )
    ) +
      ggplot2::geom_line(linewidth = 0.75) +
      ggplot2::geom_point(size = 1.5) +
      ggrepel::geom_text_repel(
        data = label_tp,
        ggplot2::aes(label = tx_label),
        size = 3,
        show.legend = FALSE,
        min.segment.length = 0.1,
        segment.alpha = 0.5,
        box.padding = 0.75,
        point.padding = 0.75,
        force = 5,
        force_pull = 0.1,
        nudge_x = 0.3,
        direction = "y",
        hjust = 0,
        max.iter = 50000,
        max.overlaps = Inf
      ) +
      ggh4x::facet_grid2(
        rows = ggplot2::vars(selected_timepoint),
        cols = ggplot2::vars(symbol),
        scales = "free_y",
        strip = ggh4x::strip_themed(
          background_y = ggh4x::elem_list_rect(
            fill = tp_color,
            color = unname(grey_palette_colors["dark"])
          ),
          text_y = ggh4x::elem_list_text(
            color = unname(grey_palette_colors["black"]),
            face = "bold",
            size = 12
          ),
          background_x = ggh4x::elem_list_rect(
            fill = tp_color,
            color = unname(grey_palette_colors["dark"])
          ),
          text_x = ggh4x::elem_list_text(
            color = unname(grey_palette_colors["black"]),
            face = "bold"
          )
        )
      ) +
      ggplot2::scale_color_manual(
        values = tx_colors,
        breaks = names(tx_colors),
        guide = "none"
      ) +
      ggplot2::scale_x_discrete(labels = c(
        C1 = "Normoxia",
        H1 = "Hypoxia (1H)",
        H3 = "Hypoxia (3H)",
        H24 = "Hypoxia (24H)"
      )) +
      ggplot2::scale_y_continuous(
        labels = scales::percent_format(accuracy = 1),
        expand = ggplot2::expansion(mult = c(0.03, 0.08))
      ) +
      ggplot2::labs(
        x = NULL,
        y = "Mean Transcript Proportion",
        title = if (show_title) title else NULL,
        subtitle = if (show_title) {
          bquote(
            "Top " * .(n_genes_per_timepoint) *
              " DTU genes per timepoint, selected by |" * Delta * " Proportion|"
          )
        } else {
          NULL
        }
      ) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        axis.title.y = ggplot2::element_text(size = 12),
        plot.title = ggplot2::element_text(face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 12)
      )
  }

  tp_order <- .drimseq_timepoint_order()

  plots <- lapply(seq_along(tp_order), function(i) {
    make_tp_plot(
      tp_raw = tp_order[[i]],
      show_title = i == 1
    )
  })

  fig_c <- patchwork::wrap_plots(
    plots,
    ncol = 1,
    heights = rel_heights
  )

  ggplot2::ggsave(
    filename = file.path(outdir, "sup_fig_3_C_drimseq_transcript_usage_trajectories.pdf"),
    plot = fig_c,
    width = 12,
    height = 10
  )

  invisible(list(
    plot = fig_c,
    by_timepoint = plots
  ))
}

# --------------------------------------------------
# Convenience Wrapper
# --------------------------------------------------

plot_sup_fig_3_all <- function(
  drimseq_results,
  outdir = NULL,
  sig_only = TRUE,
  ref_level = NULL,
  n_genes_per_timepoint = 5,
  top_tx_per_gene = 3,
  max_features_to_show = 5,
  max_effect_to_show = 1,
  label_top_n_by_features = 5,
  label_min_features = 2,
  save_tables = FALSE,
  force_recompute = FALSE,
  base_size = 15
) {
  outdir <- .set_outdir(outdir, subdir = "sup_fig_3")

  res <- run_drimseq_supplement_analysis(
    drimseq_results = drimseq_results,
    outdir = outdir,
    sig_only = sig_only,
    ref_level = ref_level,
    n_genes_per_timepoint = n_genes_per_timepoint,
    top_tx_per_gene = top_tx_per_gene,
    save_tables = save_tables,
    force_recompute = force_recompute
  )

  fig_a <- plot_drimseq_overlap_across_timepoints(
    res = res,
    outdir = outdir
  )

  fig_b <- plot_drimseq_multiplicity_vs_magnitude(
    res = res,
    outdir = outdir,
    base_size = base_size,
    max_features_to_show = max_features_to_show,
    max_effect_to_show = max_effect_to_show,
    label_top_n_by_features = label_top_n_by_features,
    label_min_features = label_min_features
  )

  fig_c <- plot_drimseq_transcript_trajectories(
    res = res,
    outdir = outdir,
    base_size = base_size
  )

  invisible(list(
    analysis = res,
    fig3_a = fig_a,
    fig3_b = fig_b,
    fig3_c = fig_c
  ))
}
