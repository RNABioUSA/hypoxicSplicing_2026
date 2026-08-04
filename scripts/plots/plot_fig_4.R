# ============================================================
# plot_fig_4.R
# ============================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
if (!nzchar(PROJECT_ROOT)) stop("PROJECT_ROOT env var not set.")

HELPERS_FILE <- file.path(PROJECT_ROOT, "scripts/utils/helpers.R")
if (!file.exists(HELPERS_FILE)) stop("helpers.R not found at: ", HELPERS_FILE)
source(HELPERS_FILE)

COLORS_FILE <- file.path(PROJECT_ROOT, "scripts/utils/color_palette.R")
if (!file.exists(COLORS_FILE)) stop("color_palette.R not found at: ", COLORS_FILE)
source(COLORS_FILE)

# -------------------------
# Internal Helpers
# -------------------------

# Return genes matching an exact TRUE/FALSE pattern across UpSet membership columns
.get_genes_matching_pattern <- function(upset_df, pattern) {
  keep <- rep(TRUE, nrow(upset_df))

  for (nm in names(pattern)) {
    keep <- keep & (upset_df[[nm]] == pattern[[nm]])
  }

  rownames(upset_df)[keep]
}

.bind_deseq_gene_full <- function(deseq_results) {
  .check_pkg("dplyr")

  if (is.null(deseq_results$results$gene_full_by_tp)) {
    stop("Expected deseq_results$results$gene_full_by_tp.")
  }

  dplyr::bind_rows(lapply(names(deseq_results$results$gene_full_by_tp), function(tp) {
    x <- as.data.frame(deseq_results$results$gene_full_by_tp[[tp]], stringsAsFactors = FALSE)
    x$timepoint <- tp
    x
  })) |>
    dplyr::select(
      timepoint,
      ensgene,
      deseq_log2fc = log2FoldChange,
      deseq_padj = padj,
      deseq_baseMean = dplyr::any_of("baseMean")
    ) |>
    dplyr::mutate(
      timepoint = factor(timepoint, levels = .timepoint_order()),
      deseq_sig = !is.na(deseq_padj) & deseq_padj < deseq_results$meta$padj_cutoff,
      abs_deseq_log2fc = abs(deseq_log2fc)
    )
}

.timepoint_order <- function() {
  c("H1", "H3", "H24")
}

.timepoint_labels <- function() {
  stats::setNames(
    names(timepoint_base_colors)[grepl("Hypoxia", names(timepoint_base_colors))],
    .timepoint_order()
  )
}

.splicing_vs_deseq_colors <- function(tool_name) {
  tool_base <- splicing_tool_colors[tool_name]

  if (is.na(tool_base) || length(tool_base) == 0) {
    stop("No splicing_tool_colors entry found for tool_name: ", tool_name)
  }

  tool_shades <- .generate_shaded_palette(
    base_colors = tool_base,
    labels = c("light", "base", "dark")
  )

  c(
    "Splicing Level Dominant" = unname(tool_shades[[tool_name]][["dark"]]),
    "Both Levels Dominant" = unname(grey_palette_colors[["dark"]]),
    "Gene Level Dominant" = unname(grey_palette_colors[["mid"]]),
    "Low Effect" = unname(tool_shades[[tool_name]][["base"]])
  )
}

.splicing_vs_deseq_categories <- function(
  df,
  y_col,
  lfc_cutoff,
  effect_cutoff
) {
  .check_pkg("dplyr")

  category_levels <- c(
    "Splicing Level Dominant",
    "Both Levels Dominant",
    "Gene Level Dominant",
    "Low Effect"
  )

  df |>
    dplyr::mutate(
      category = dplyr::case_when(
        abs_deseq_log2fc >= lfc_cutoff & .data[[y_col]] >= effect_cutoff ~
          "Both Levels Dominant",
        abs_deseq_log2fc < lfc_cutoff & .data[[y_col]] >= effect_cutoff ~
          "Splicing Level Dominant",
        abs_deseq_log2fc >= lfc_cutoff & .data[[y_col]] < effect_cutoff ~
          "Gene Level Dominant",
        TRUE ~ "Low Effect"
      ),
      category = factor(category, levels = category_levels)
    )
}

# --------------------------------------------------
# Splicing Tool Data Builders
# --------------------------------------------------

make_dexseq_vs_deseq_df <- function(
  dexseq_results,
  deseq_results,
  padj_cutoff = 0.10
) {
  .check_pkg(c("dplyr"))

  if (!exists(".make_dexseq_multiplicity_magnitude_df", mode = "function")) {
    stop(
      "Expected helper function .make_dexseq_multiplicity_magnitude_df() ",
      "to be available. Source the DEXSeq figure script first."
    )
  }

  dex_gene <- .make_dexseq_multiplicity_magnitude_df(
    dexseq_results = dexseq_results,
    padj_cutoff = padj_cutoff,
    sig_only = TRUE
  ) |>
    dplyr::rename(
      n_sig_features = n_sig_exons,
      max_abs_effect = max_abs_exon_effect,
      mean_abs_effect = mean_abs_exon_effect
    )

  deseq_df <- .bind_deseq_gene_full(deseq_results)

  dex_gene |>
    dplyr::left_join(
      deseq_df,
      by = c("timepoint", "ensgene")
    ) |>
    dplyr::mutate(
      timepoint = factor(timepoint, levels = .timepoint_order()),
      symbol = dplyr::if_else(is.na(symbol) | symbol == "", ensgene, symbol),
      abs_deseq_log2fc = abs(deseq_log2fc)
    )
}

make_drimseq_vs_deseq_df <- function(
  drimseq_results,
  deseq_results,
  ref_level = NULL
) {
  .check_pkg(c("dplyr"))

  if (!exists("summarize_drimseq_delta_proportions", mode = "function")) {
    stop("Expected helper function summarize_drimseq_delta_proportions() to be available.")
  }

  if (is.null(ref_level)) {
    ref_level <- drimseq_results$meta$ref_level
  }

  dprop_df <- summarize_drimseq_delta_proportions(
    drimseq_results,
    ref_level = ref_level
  )

  drim_gene_df <- dprop_df |>
    dplyr::filter(significant %in% TRUE, is.finite(abs_dprop)) |>
    dplyr::group_by(timepoint, ensgene, symbol) |>
    dplyr::summarize(
      n_sig_features = dplyr::n_distinct(feature_id),
      max_abs_effect = max(abs_dprop, na.rm = TRUE),
      mean_abs_effect = mean(abs_dprop, na.rm = TRUE),
      min_tool_padj = suppressWarnings(min(padj, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      min_tool_padj = dplyr::if_else(is.infinite(min_tool_padj), NA_real_, min_tool_padj)
    )

  deseq_df <- .bind_deseq_gene_full(deseq_results)

  drim_gene_df |>
    dplyr::left_join(deseq_df, by = c("timepoint", "ensgene")) |>
    dplyr::mutate(
      timepoint = factor(timepoint, levels = .timepoint_order()),
      symbol = dplyr::if_else(is.na(symbol) | symbol == "", ensgene, symbol),
      abs_deseq_log2fc = abs(deseq_log2fc)
    )
}

make_suppa_vs_deseq_df <- function(suppa_results, deseq_results) {
  .check_pkg(c("dplyr"))

  if (is.null(suppa_results$results$event_sig_all)) {
    stop("Expected suppa_results$results$event_sig_all.")
  }

  suppa_gene <- suppa_results$results$event_sig_all |>
    dplyr::filter(!is.na(pvalue), is.finite(dpsi)) |>
    dplyr::group_by(timepoint, ensgene, symbol) |>
    dplyr::summarize(
      n_sig_features = dplyr::n_distinct(event_id),
      max_abs_effect = max(abs(dpsi), na.rm = TRUE),
      mean_abs_effect = mean(abs(dpsi), na.rm = TRUE),
      min_tool_pvalue = min(pvalue, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      min_tool_pvalue = dplyr::if_else(is.infinite(min_tool_pvalue), NA_real_, min_tool_pvalue)
    )

  deseq_df <- .bind_deseq_gene_full(deseq_results)

  suppa_gene |>
    dplyr::left_join(deseq_df, by = c("timepoint", "ensgene")) |>
    dplyr::mutate(
      timepoint = factor(timepoint, levels = .timepoint_order()),
      symbol = dplyr::if_else(is.na(symbol) | symbol == "", ensgene, symbol),
      abs_deseq_log2fc = abs(deseq_log2fc)
    )
}

# =========================================================
# Main Analysis Functions
# =========================================================

run_gene_vs_transcript_expression_analysis <- function(
  deseq_results,
  dexseq_results = NULL,
  drimseq_results = NULL,
  suppa_results = NULL,
  drimseq_ref_level = NULL
) {
  .check_pkg("dplyr")

  out <- list()

  if (!is.null(dexseq_results)) {
    dex_df <- make_dexseq_vs_deseq_df(
      dexseq_results = dexseq_results,
      deseq_results = deseq_results
    ) |>
      dplyr::rename(
        n_sig_features = dplyr::any_of("n_sig_exons"),
        max_abs_effect = dplyr::any_of("max_abs_effect"),
        mean_abs_effect = dplyr::any_of("mean_abs_effect")
      ) |>
      dplyr::mutate(tool = "DEXSeq")

    out$dexseq <- dex_df
  }

  if (!is.null(drimseq_results)) {
    out$drimseq <- make_drimseq_vs_deseq_df(
      drimseq_results = drimseq_results,
      deseq_results = deseq_results,
      ref_level = drimseq_ref_level
    ) |>
      dplyr::mutate(tool = "DRIMSeq")
  }

  if (!is.null(suppa_results)) {
    out$suppa2 <- make_suppa_vs_deseq_df(
      suppa_results = suppa_results,
      deseq_results = deseq_results
    ) |>
      dplyr::mutate(tool = "SUPPA2")
  }

  if (length(out) == 0) {
    stop("No splicing result objects were provided.")
  }

  invisible(list(
    by_tool = out,
    combined = dplyr::bind_rows(out)
  ))
}

run_gene_vs_transcript_overlap_analysis <- function(dexseq_results,
                                                    drimseq_results,
                                                    suppa_results,
                                                    deseq_results,
                                                    timepoints = c("H1", "H3", "H24")) {
  .check_pkg("ComplexUpset")
  .check_pkg("ComplexHeatmap")

  deg_genes_by_tp <- .get_sig_genes_by_tp(deseq_results$results$gene_sig_by_tp)
  deg_genes_by_tp <- .subset_gene_lists_by_tp(deg_genes_by_tp, timepoints)

  dex_genes_by_tp <- .get_sig_genes_by_tp(dexseq_results$results$gene_sig_by_tp)
  drim_genes_by_tp <- .get_sig_genes_by_tp(drimseq_results$results$gene_sig_by_tp)
  suppa_genes_by_tp <- .get_sig_genes_by_tp(suppa_results$results$gene_sig_by_tp)

  dex_genes_by_tp <- .subset_gene_lists_by_tp(dex_genes_by_tp, timepoints)
  drim_genes_by_tp <- .subset_gene_lists_by_tp(drim_genes_by_tp, timepoints)
  suppa_genes_by_tp <- .subset_gene_lists_by_tp(suppa_genes_by_tp, timepoints)

  as_genes_by_tp <- .union_gene_lists_by_tp(
    dex_genes_by_tp,
    drim_genes_by_tp,
    suppa_genes_by_tp
  )
  as_genes_by_tp <- .subset_gene_lists_by_tp(as_genes_by_tp, timepoints)

  upset_list <- c(
    stats::setNames(deg_genes_by_tp, paste0(timepoints, "_DEG")),
    stats::setNames(as_genes_by_tp, paste0(timepoints, "_ASG"))
  )

  upset_df <- as.data.frame(list_to_matrix(upset_list))
  groups <- colnames(upset_df)
  upset_df[groups] <- upset_df[groups] == 1

  # Exact highlighted patterns used for both UpSet and alluvial
  highlight_patterns <- list(
    "Group 1" = c(H1_DEG = FALSE, H1_ASG = TRUE, H3_DEG = FALSE, H3_ASG = TRUE, H24_DEG = TRUE, H24_ASG = FALSE),
    "Group 2" = c(H1_DEG = FALSE, H1_ASG = TRUE, H3_DEG = TRUE, H3_ASG = FALSE, H24_DEG = TRUE, H24_ASG = FALSE),
    "Group 3" = c(H1_DEG = FALSE, H1_ASG = TRUE, H3_DEG = TRUE, H3_ASG = FALSE, H24_DEG = FALSE, H24_ASG = TRUE),
    "Group 4" = c(H1_DEG = FALSE, H1_ASG = TRUE, H3_DEG = FALSE, H3_ASG = TRUE, H24_DEG = FALSE, H24_ASG = FALSE),
    "Group 5" = c(H1_DEG = FALSE, H1_ASG = TRUE, H3_DEG = FALSE, H3_ASG = FALSE, H24_DEG = FALSE, H24_ASG = TRUE),
    "Group 6" = c(H1_DEG = FALSE, H1_ASG = TRUE, H3_DEG = FALSE, H3_ASG = FALSE, H24_DEG = TRUE, H24_ASG = FALSE)
  )

  highlighted_sets <- lapply(highlight_patterns, function(pattern) {
    .get_genes_matching_pattern(upset_df, pattern)
  })

  alluvial_df <- data.frame(
    H1 = rep("H1_ASG", length(highlight_patterns)),
    H3 = c("H3_ASG", "H3_DEG", "H3_DEG", "H3_ASG", "H3_NSG", "H3_NSG"),
    H24 = c("H24_DEG", "H24_DEG", "H24_ASG", "H24_NSG", "H24_ASG", "H24_DEG"),
    Group = names(highlight_patterns),
    Freq = vapply(highlighted_sets, length, integer(1)),
    stringsAsFactors = FALSE
  )

  out <- list(
    inputs = list(
      timepoints = timepoints
    ),
    results = list(
      deg_genes_by_tp = deg_genes_by_tp,
      as_genes_by_tp = as_genes_by_tp,
      dex_genes_by_tp = dex_genes_by_tp,
      drim_genes_by_tp = drim_genes_by_tp,
      suppa_genes_by_tp = suppa_genes_by_tp,
      upset_list = upset_list,
      upset_df = upset_df,
      groups = groups,
      highlight_patterns = highlight_patterns,
      highlighted_sets = highlighted_sets,
      alluvial_df = alluvial_df
    )
  )

  out
}

# =========================================================
# Figure 4A - Scatter Plot
# =========================================================

plot_splicing_vs_deseq_scatter <- function(
  df,
  y_col,
  tool_name,
  x_limits = NULL,
  y_limits = NULL,
  y_breaks = NULL,
  y_label,
  lfc_cutoff = 1,
  effect_cutoff,
  label_top_n = 5,
  base_size = 15,
  title = NULL
) {
  .check_pkg(c("dplyr", "ggplot2", "ggrepel", "ggh4x"))

  plot_df <- df |>
    dplyr::filter(is.finite(deseq_log2fc), is.finite(.data[[y_col]])) |>
    .splicing_vs_deseq_categories(
      y_col = y_col,
      lfc_cutoff = lfc_cutoff,
      effect_cutoff = effect_cutoff
    )

  label_df <- plot_df |>
    dplyr::group_by(timepoint) |>
    dplyr::arrange(dplyr::desc(.data[[y_col]]), abs_deseq_log2fc, .by_group = TRUE) |>
    dplyr::slice_head(n = label_top_n) |>
    dplyr::ungroup()

  plot_df$tool <- tool_name
  label_df$tool <- tool_name

  if (is.null(x_limits)) {
    max_abs_x <- max(abs(plot_df$deseq_log2fc), na.rm = TRUE)
    max_abs_x <- max(max_abs_x, lfc_cutoff)
    x_limits <- ceiling(max_abs_x)
  }

  p_scatter <- ggplot2::ggplot(plot_df, ggplot2::aes(x = deseq_log2fc, y = .data[[y_col]])) +
    ggplot2::geom_vline(
      xintercept = c(-lfc_cutoff, lfc_cutoff),
      linetype = "dashed",
      linewidth = 0.4,
      color = unname(grey_palette_colors["dark"])
    ) +
    ggplot2::geom_hline(
      yintercept = effect_cutoff,
      linetype = "dashed",
      linewidth = 0.4,
      color = unname(grey_palette_colors["dark"])
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = category),
      alpha = 0.75,
      size = 1.75
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(label = symbol),
      size = 3,
      max.overlaps = Inf,
      max.iter = 50000,
      max.time = 5,
      box.padding = 0.8,
      point.padding = 0.8,
      min.segment.length = 0.1,
      segment.alpha = 0.5,
      show.legend = FALSE
    ) +
    ggh4x::facet_grid2(
      tool ~ timepoint,
      labeller = ggplot2::labeller(
        timepoint = .timepoint_labels()
      ),
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = unname(timepoint_base_colors[.timepoint_labels()]),
          color = unname(grey_palette_colors["black"])
        ),
        text_x = ggh4x::elem_list_text(
          color = "white",
          face = "bold"
        ),
        background_y = ggh4x::elem_list_rect(
          fill = unname(splicing_tool_colors[[tool_name]]),
          color = unname(grey_palette_colors["black"])
        ),
        text_y = ggh4x::elem_list_text(
          color = "white",
          face = "bold"
        )
      )
    ) +
    ggplot2::scale_color_manual(
      values = .splicing_vs_deseq_colors(tool_name),
      breaks = names(.splicing_vs_deseq_colors(tool_name)),
      drop = FALSE
    ) +
    ggplot2::labs(
      x = expression("Gene-Level log"[2] * "FC"),
      y = y_label,
      color = NULL,
      title = title
    ) +
    ggplot2::coord_cartesian(xlim = c(-x_limits, x_limits)) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(y_limits)) {
    p_scatter <- p_scatter +
      ggplot2::scale_y_continuous(
        limits = y_limits,
        breaks = y_breaks
      )
  }

  p_scatter
}

plot_tool_row_title <- function(
  tool_name,
  base_size = 18,
  fill = "white"
) {
  .check_pkg("ggplot2")

  tool_col <- .tool_color(tool_name)

  ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 0.5,
      label = tool_name,
      hjust = 0,
      fontface = "bold",
      color = tool_col,
      size = base_size / 3
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0, 1),
      clip = "off"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(
        fill = fill,
        color = NA
      ),
      plot.margin = ggplot2::margin(2, 2, 0, 2)
    )
}

plot_fig_4_A_tool_panel <- function(
  plot_data,
  tool_key,
  tool_name,
  y_label,
  effect_cutoff,
  x_limit = 3,
  y_limit = NULL,
  y_breaks = NULL,
  lfc_cutoff = 1,
  label_top_n = 5,
  base_size = 15
) {
  if (is.null(plot_data$by_tool[[tool_key]])) {
    stop("plot_data$by_tool$", tool_key, " not found.")
  }

  plot_splicing_vs_deseq_scatter(
    df = plot_data$by_tool[[tool_key]],
    y_col = "max_abs_effect",
    tool_name = tool_name,
    y_label = y_label,
    x_limits = x_limit,
    y_limits = y_limit,
    y_breaks = y_breaks,
    lfc_cutoff = lfc_cutoff,
    effect_cutoff = effect_cutoff,
    label_top_n = label_top_n,
    base_size = base_size,
    title = NULL
  )
}

plot_fig_4_A <- function(
  scatter_data,
  outdir = NULL,
  filename = "fig_4_A_gene_vs_transcript_scatter.pdf",
  width = 12,
  height = 18,
  base_size = 15,
  save_pdf = TRUE
) {
  .check_pkg(c("patchwork", "ggplot2"))

  outdir <- .set_outdir(outdir, subdir = "fig_4")

  p_dex <- plot_fig_4_A_tool_panel(
    plot_data = scatter_data,
    tool_key = "dexseq",
    tool_name = "DEXSeq",
    y_label = expression("Max |Exon log"[2] * "FC| per Gene"),
    effect_cutoff = 1,
    x_limit = 3,
    base_size = base_size
  )

  p_drim <- plot_fig_4_A_tool_panel(
    plot_data = scatter_data,
    tool_key = "drimseq",
    tool_name = "DRIMSeq",
    y_label = expression("Max |" * Delta * " Transcript Proportion| per Gene"),
    effect_cutoff = 0.20,
    x_limit = 3,
    y_limit = c(0, 1),
    y_breaks = seq(0, 1, 0.2),
    base_size = base_size
  )

  p_suppa <- plot_fig_4_A_tool_panel(
    plot_data = scatter_data,
    tool_key = "suppa2",
    tool_name = "SUPPA2",
    y_label = expression("Max |" * Delta * "PSI| per Gene"),
    effect_cutoff = 0.20,
    x_limit = 3,
    y_limit = c(0, 1),
    y_breaks = seq(0, 1, 0.2),
    base_size = base_size
  )

  dex_row <- (
    patchwork::wrap_elements(full = plot_tool_row_title("DEXSeq")) /
      patchwork::wrap_elements(full = p_dex)
  ) +
    patchwork::plot_layout(heights = c(0.08, 1))

  drim_row <- (
    patchwork::wrap_elements(full = plot_tool_row_title("DRIMSeq")) /
      patchwork::wrap_elements(full = p_drim)
  ) +
    patchwork::plot_layout(heights = c(0.08, 1))

  suppa_row <- (
    patchwork::wrap_elements(full = plot_tool_row_title("SUPPA2")) /
      patchwork::wrap_elements(full = p_suppa)
  ) +
    patchwork::plot_layout(heights = c(0.08, 1))

  fig4_a <- patchwork::wrap_plots(
    dex_row,
    drim_row,
    suppa_row,
    ncol = 1,
    heights = c(1, 1, 1)
  )

  if (isTRUE(save_pdf)) {
    ggplot2::ggsave(
      filename = file.path(outdir, filename),
      plot = fig4_a,
      width = width,
      height = height
    )
  }

  invisible(fig4_a)
}

# =========================================================
# Figure 4B - UpSet plot
# =========================================================

plot_gene_vs_transcript_upset <- function(overlap_data,
                                          min_size = 10,
                                          base_size = 15) {
  .check_pkg(c("ComplexUpset", "ggplot2", "scales"))

  timepoints <- overlap_data$inputs$timepoints
  upset_df <- overlap_data$results$upset_df
  groups <- overlap_data$results$groups
  highlight_patterns <- overlap_data$results$highlight_patterns

  set_labels <- c(
    H1_DEG = "1H - DEG",
    H3_DEG = "3H  - DEG",
    H24_DEG = "24H - DEG",
    H1_ASG = "1H - ASG",
    H3_ASG = "3H - ASG",
    H24_ASG = "24H - ASG"
  )

  pretty_labeller <- function(sets) set_labels[sets]

  tp_fill_colors <- .get_timepoint_shaded_colors(
    timepoints = timepoints,
    labels = c("ASG", "OVERLAP", "DEG")
  )

  set_colors <- c(
    H1_DEG = unname(tp_fill_colors$H1["DEG"]),
    H3_DEG = unname(tp_fill_colors$H3["DEG"]),
    H24_DEG = unname(tp_fill_colors$H24["DEG"]),
    H1_ASG = unname(tp_fill_colors$H1["ASG"]),
    H3_ASG = unname(tp_fill_colors$H3["ASG"]),
    H24_ASG = unname(tp_fill_colors$H24["ASG"])
  )

  set_queries <- list(
    ComplexUpset::upset_query(set = "H1_DEG", fill = set_colors["H1_DEG"]),
    ComplexUpset::upset_query(set = "H1_ASG", fill = set_colors["H1_ASG"]),
    ComplexUpset::upset_query(set = "H24_ASG", fill = set_colors["H24_ASG"]),
    ComplexUpset::upset_query(set = "H3_DEG", fill = set_colors["H3_DEG"]),
    ComplexUpset::upset_query(set = "H24_DEG", fill = set_colors["H24_DEG"]),
    ComplexUpset::upset_query(set = "H3_ASG", fill = set_colors["H3_ASG"])
  )

  group_colors <- c(
    "Group 1" = unname(category_base_colors["Category 1"]),
    "Group 2" = unname(category_base_colors["Category 4"]),
    "Group 3" = unname(category_base_colors["Category 3"]),
    "Group 4" = unname(category_base_colors["Category 2"]),
    "Group 5" = unname(category_base_colors["Category 5"]),
    "Group 6" = unname(category_base_colors["Category 6"])
  )

  intersection_queries <- lapply(names(highlight_patterns), function(group_name) {
    pattern <- highlight_patterns[[group_name]]
    active_sets <- names(pattern)[pattern]

    ComplexUpset::upset_query(
      intersect = active_sets,
      color = group_colors[[group_name]],
      fill = group_colors[[group_name]],
      only_components = c("intersections_matrix", "Significant Genes")
    )
  })

  set_metadata <- data.frame(
    set = groups,
    set_id = groups,
    stringsAsFactors = FALSE
  )

  set_metadata$set_id <- factor(set_metadata$set_id, levels = groups)

  stripe_cols <- vapply(
    set_colors[groups],
    function(x) .blend_color(x, "white", 0.60),
    character(1)
  )
  names(stripe_cols) <- groups

  p_upset <- ComplexUpset::upset(
    upset_df,
    groups,
    min_size = min_size,
    name = "Gene Set Overlap",
    base_annotations = list(
      "Significant Genes" =
        ComplexUpset::intersection_size(
          mapping = ggplot2::aes(fill = "bars_color"),
          text = list(vjust = -0.5, hjust = -0.05, size = 3, angle = 30),
          counts = TRUE,
          width = 0.75
        ) +
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
          ggplot2::scale_fill_manual(
            guide = "none",
            values = c("bars_color" = unname(grey_palette_colors["black"]))
          )
    ),
    themes = ComplexUpset::upset_default_themes(
      text = ggplot2::element_text(size = base_size)
    ),
    stripes = ComplexUpset::upset_stripes(
      mapping = ggplot2::aes(color = set_id),
      data = set_metadata,
      colors = stripe_cols
    ),
    labeller = pretty_labeller,
    set_sizes = ComplexUpset::upset_set_size(
      geom = ggplot2::geom_bar(ggplot2::aes(fill = "bars_color"))
    ) +
      ggplot2::scale_fill_manual(
        values = c("bars_color" = unname(grey_palette_colors["black"])),
        guide = "none"
      ) +
      ggplot2::geom_text(
        stat = "count",
        ggplot2::aes(label = ggplot2::after_stat(count)),
        vjust = 0.5,
        hjust = -0.1,
        size = 3
      ) +
      ggplot2::ylab("Gene Set Size"),
    queries = c(set_queries, intersection_queries),
    matrix = ComplexUpset::intersection_matrix(
      geom = ggplot2::geom_point(size = 3),
      segment = ggplot2::geom_segment(linewidth = 1.5),
      outline_color = list(
        active = "transparent",
        inactive = unname(grey_palette_colors["light"])
      )
    )
  ) +
    ggplot2::labs(
      title = "Hypoxia Gene Sets"
    ) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 15, face = "bold")
    )

  p_upset
}

plot_fig_4_B <- function(
  overlap_data,
  outdir = NULL,
  filename = "fig_4_B_gene_vs_transcript_upset.pdf",
  min_size = 10,
  width = 15,
  height = 10,
  base_size = 15,
  save_pdf = TRUE
) {
  outdir <- .set_outdir(outdir, subdir = "fig_4")

  fig4_b <- plot_gene_vs_transcript_upset(
    overlap_data = overlap_data,
    min_size = min_size,
    base_size = base_size
  )

  if (isTRUE(save_pdf)) {
    ggplot2::ggsave(
      filename = file.path(outdir, filename),
      plot = fig4_b,
      width = width,
      height = height
    )
  }

  invisible(fig4_b)
}

# =========================================================
# Figure 4C - Alluvial Plot
# =========================================================

plot_gene_vs_transcript_alluvial <- function(overlap_data,
                                             base_size = 15) {
  .check_pkg(c("ggplot2", "ggalluvial", "ggnewscale"))

  df <- overlap_data$results$alluvial_df

  ribbon_order <- c("Group 1", "Group 4", "Group 3", "Group 2", "Group 5", "Group 6")
  df$Group <- factor(df$Group, levels = ribbon_order)

  alluvial_label_df <- df %>%
    dplyr::arrange(dplyr::desc(Group)) %>%
    dplyr::mutate(
      percent = Freq / sum(Freq),
      percent_label = scales::percent(percent, accuracy = 0.1),
      ymax = cumsum(Freq),
      ymin = ymax - Freq,
      y = (ymin + ymax) / 2,
      x = 1.5
    )

  tp_fill_colors <- .get_timepoint_shaded_colors(
    timepoints = c("H1", "H3", "H24"),
    labels = c("ASG", "NSG", "DEG")
  )

  stratum_colors <- c(
    H1_ASG = unname(tp_fill_colors$H1["ASG"]),
    H3_ASG = unname(tp_fill_colors$H3["ASG"]),
    H3_DEG = unname(tp_fill_colors$H3["DEG"]),
    H3_NSG = unname(tp_fill_colors$H3["NSG"]),
    H24_ASG = unname(tp_fill_colors$H24["ASG"]),
    H24_DEG = unname(tp_fill_colors$H24["DEG"]),
    H24_NSG = unname(tp_fill_colors$H24["NSG"])
  )

  group_colors <- c(
    "Group 1" = unname(category_base_colors["Category 1"]),
    "Group 2" = unname(category_base_colors["Category 4"]),
    "Group 3" = unname(category_base_colors["Category 3"]),
    "Group 4" = unname(category_base_colors["Category 2"]),
    "Group 5" = unname(category_base_colors["Category 5"]),
    "Group 6" = unname(category_base_colors["Category 6"])
  )

  df_long <- ggalluvial::to_lodes_form(
    df,
    key = "Timepoint",
    axes = 1:3
  )

  stratum_labels <- c(
    H1_ASG = "ASG",
    H3_ASG = "ASG",
    H3_DEG = "DEG",
    H3_NSG = "NSG",
    H24_ASG = "ASG",
    H24_DEG = "DEG",
    H24_NSG = "NSG"
  )

  p_alluvial <- ggplot2::ggplot(
    data = df_long,
    ggplot2::aes(
      x = Timepoint,
      stratum = stratum,
      alluvium = alluvium,
      y = Freq
    )
  ) +
    ggplot2::scale_x_discrete(
      limits = c("H1", "H3", "H24"),
      labels = c(
        H1 = "1H",
        H3 = "3H",
        H24 = "24H"
      ),
      expand = c(0.12, 0.05)
    ) +
    ggalluvial::geom_alluvium(
      ggplot2::aes(fill = Group),
      alpha = 0.85,
      width = 0.20,
      knot.pos = 0.5,
      color = "black",
      linewidth = 0.3
    ) +
    ggplot2::scale_fill_manual(
      values = group_colors,
      guide = "none"
    ) +
    ggnewscale::new_scale_fill() +
    ggalluvial::geom_stratum(
      ggplot2::aes(fill = stratum),
      width = 0.30,
      color = "black",
      linewidth = 0.3
    ) +
    ggplot2::scale_fill_manual(
      values = stratum_colors,
      guide = "none"
    ) +
    ggalluvial::stat_stratum(
      ggplot2::aes(label = ggplot2::after_stat(stratum_labels[stratum])),
      geom = "text",
      size = 3.5,
      lineheight = 0.95
    ) +
    ggplot2::labs(
      x = "Hypoxia",
      y = "Gene Count"
    ) +
    ggplot2::geom_text(
      data = alluvial_label_df,
      ggplot2::aes(x = x, y = y, label = percent_label),
      inherit.aes = FALSE,
      size = 3.5
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none"
    )

  p_alluvial
}

plot_fig_4_C <- function(
  overlap_data,
  outdir = NULL,
  filename = "fig_4_C_gene_vs_transcript_alluvial.pdf",
  width = 9,
  height = 6,
  base_size = 15,
  save_pdf = TRUE
) {
  outdir <- .set_outdir(outdir, subdir = "fig_4")

  fig4_c <- plot_gene_vs_transcript_alluvial(
    overlap_data = overlap_data,
    base_size = base_size
  )

  if (isTRUE(save_pdf)) {
    ggplot2::ggsave(
      filename = file.path(outdir, filename),
      plot = fig4_c,
      width = width,
      height = height
    )
  }

  invisible(fig4_c)
}

# =========================================================
# Convenience Wrapper
# =========================================================

run_fig_4_analysis <- function(
  deseq_results,
  dexseq_results,
  drimseq_results,
  suppa_results,
  outdir = NULL,
  force = FALSE
) {
  outdir <- .set_outdir(outdir, subdir = "fig_4")

  scatter_rds <- file.path(outdir, "fig_4_gene_vs_transcript_expression_analysis.rds")
  overlap_rds <- file.path(outdir, "fig_4_gene_vs_transcript_overlap_analysis.rds")

  if (file.exists(scatter_rds) && !force) {
    message("Loading cached expression data: ", scatter_rds)
    scatter_data <- readRDS(scatter_rds)
  } else {
    scatter_data <- run_gene_vs_transcript_expression_analysis(
      deseq_results = deseq_results,
      dexseq_results = dexseq_results,
      drimseq_results = drimseq_results,
      suppa_results = suppa_results
    )
    saveRDS(scatter_data, scatter_rds)
  }

  if (file.exists(overlap_rds) && !force) {
    message("Loading cached overlap data: ", overlap_rds)
    overlap_data <- readRDS(overlap_rds)
  } else {
    overlap_data <- run_gene_vs_transcript_overlap_analysis(
      dexseq_results = dexseq_results,
      drimseq_results = drimseq_results,
      suppa_results = suppa_results,
      deseq_results = deseq_results
    )
    saveRDS(overlap_data, overlap_rds)
  }

  invisible(list(
    scatter = scatter_data,
    overlap = overlap_data,
    paths = list(
      outdir = outdir,
      scatter_rds = scatter_rds,
      overlap_rds = overlap_rds
    )
  ))
}

plot_fig_4_all <- function(
  deseq_results,
  dexseq_results,
  drimseq_results,
  suppa_results,
  outdir = NULL,
  force_analysis = FALSE,
  save_pdf = TRUE
) {
  outdir <- .set_outdir(outdir, subdir = "fig_4")

  analysis <- run_fig_4_analysis(
    deseq_results = deseq_results,
    dexseq_results = dexseq_results,
    drimseq_results = drimseq_results,
    suppa_results = suppa_results,
    outdir = outdir,
    force = force_analysis
  )

  fig_a <- plot_fig_4_A(
    scatter_data = analysis$scatter,
    outdir = outdir,
    save_pdf = save_pdf
  )

  fig_b <- plot_fig_4_B(
    overlap_data = analysis$overlap,
    outdir = outdir,
    save_pdf = save_pdf
  )

  fig_c <- plot_fig_4_C(
    overlap_data = analysis$overlap,
    outdir = outdir,
    save_pdf = save_pdf
  )

  invisible(list(
    analysis = analysis,
    plots = list(
      A = fig_a,
      B = fig_b,
      C = fig_c
    ),
    paths = list(outdir = outdir)
  ))
}
