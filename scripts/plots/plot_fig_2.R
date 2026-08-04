# ============================================================
# plot_fig_2.R
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
# Transcript-Level Analysis
# --------------------------------------------------

run_transcript_level_analysis <- function(
  deseq_results,
  dexseq_results,
  drimseq_results,
  suppa_results,
  outdir = NULL,
  save_tables = FALSE
) {
  .check_pkg("dplyr")

  if (is.null(outdir)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_2")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  outfile <- file.path(outdir, "fig_2_splicing_vs_deg_analysis.rds")

  if (file.exists(outfile)) {
    message("Loading cached transcript-level analysis from: ", outfile)
    return(readRDS(outfile))
  }

  timepoints <- c("H1", "H3", "H24")
  tool_order <- c("DEXSeq", "DRIMSeq", "SUPPA2")

  deg_genes_by_tp <- .subset_gene_lists_by_tp(
    .get_sig_genes_by_tp(deseq_results$results$gene_sig_by_tp),
    timepoints = timepoints
  )

  dex_genes_by_tp <- .subset_gene_lists_by_tp(
    .get_sig_genes_by_tp(dexseq_results$results$gene_sig_by_tp),
    timepoints = timepoints
  )

  drim_genes_by_tp <- .subset_gene_lists_by_tp(
    .get_sig_genes_by_tp(drimseq_results$results$gene_sig_by_tp),
    timepoints = timepoints
  )

  suppa_genes_by_tp <- .subset_gene_lists_by_tp(
    .get_sig_genes_by_tp(suppa_results$results$gene_sig_by_tp),
    timepoints = timepoints
  )

  splicing_genes_by_tp <- .union_gene_lists_by_tp(
    dex_genes_by_tp,
    drim_genes_by_tp,
    suppa_genes_by_tp
  )
  splicing_genes_by_tp <- .subset_gene_lists_by_tp(splicing_genes_by_tp, timepoints = timepoints)

  # fig2_c data
  overlap_bar_df <- do.call(rbind, lapply(timepoints, function(tp) {
    deg_genes <- unique(deg_genes_by_tp[[tp]])

    tool_gene_sets <- list(
      DEXSeq = unique(dex_genes_by_tp[[tp]]),
      DRIMSeq = unique(drim_genes_by_tp[[tp]]),
      SUPPA2 = unique(suppa_genes_by_tp[[tp]])
    )

    do.call(rbind, lapply(names(tool_gene_sets), function(tool) {
      asg_genes <- unique(tool_gene_sets[[tool]])
      overlap_genes <- intersect(deg_genes, asg_genes)

      n_deg <- length(deg_genes)
      n_asg <- length(asg_genes)
      n_overlap <- length(overlap_genes)

      data.frame(
        timepoint = tp,
        method = tool,
        metric = c(
          "% DEGs Also Identified as Alternatively Spliced",
          "% ASGs Also Identified as Differentially Expressed"
        ),
        value = c(
          if (n_deg > 0) n_overlap / n_deg else NA_real_,
          if (n_asg > 0) n_overlap / n_asg else NA_real_
        ),
        stringsAsFactors = FALSE
      )
    }))
  }))

  # fig2_d data
  overlap_compare_df <- do.call(rbind, lapply(timepoints, function(tp) {
    deg_genes <- unique(deg_genes_by_tp[[tp]])
    asg_genes <- unique(splicing_genes_by_tp[[tp]])

    shared_genes <- intersect(deg_genes, asg_genes)
    deg_only_genes <- setdiff(deg_genes, asg_genes)
    asg_only_genes <- setdiff(asg_genes, deg_genes)

    data.frame(
      timepoint = tp,
      group = c("ASG Genes", "Shared Genes", "DEG Genes"),
      count = c(length(asg_only_genes), length(shared_genes), length(deg_only_genes)),
      stringsAsFactors = FALSE
    )
  }))

  out <- list(
    meta = list(
      outdir = outdir,
      timepoints = timepoints,
      tool_order = tool_order
    ),
    gene_sets = list(
      deg_genes_by_tp = deg_genes_by_tp,
      dex_genes_by_tp = dex_genes_by_tp,
      drim_genes_by_tp = drim_genes_by_tp,
      suppa_genes_by_tp = suppa_genes_by_tp,
      splicing_genes_by_tp = splicing_genes_by_tp
    ),
    overlap = list(
      overlap_bar_df = overlap_bar_df,
      compare_df = overlap_compare_df
    )
  )

  if (isTRUE(save_tables)) {
    utils::write.table(
      overlap_bar_df,
      file = file.path(outdir, "fig_2_C_deg_splicing_overlap.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    utils::write.table(
      overlap_compare_df,
      file = file.path(outdir, "fig_2_D_deg_vs_splicing_gene_overlap.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  saveRDS(out, file.path(outdir, "fig_2_splicing_vs_deg_analysis.rds"))
  out
}

# --------------------------------------------------
# Figure 2A - Splicing Tool Example Plots
# --------------------------------------------------

plot_splicing_tool_examples <- function(
  base_size = 15,
  outdir = NULL,
  title = NULL
) {
  .check_pkg("ggplot2")

  if (is.null(outdir)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_2")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  example_colors <- .generate_shaded_palette(
    base_colors = splicing_tool_colors["Other"],
    labels = c("light", "base", "dark")
  )

  condition_colors <- c(
    Normoxia = unname(grey_palette_colors["black"]),
    Hypoxia  = unname(example_colors$Other["base"])
  )

  base_theme <- ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.position = "right"
    )

  # DEXSeq example
  dex_data <- data.frame(
    Exon = rep(c("Exon 1", "Exon 2"), each = 2),
    Condition = rep(c("Normoxia", "Hypoxia"), times = 2),
    Count = c(60, 70, 65, 10),
    stringsAsFactors = FALSE
  )
  dex_data$Condition <- factor(dex_data$Condition, levels = c("Normoxia", "Hypoxia"))

  plot_dex <- ggplot2::ggplot(
    dex_data,
    ggplot2::aes(x = Exon, y = Count, fill = Condition)
  ) +
    ggplot2::geom_col(position = "dodge", width = 0.75) +
    ggplot2::scale_fill_manual(values = condition_colors) +
    ggplot2::scale_y_continuous(
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100)
    ) +
    ggplot2::labs(
      x = "",
      y = "Exon Counts",
      title = "DEXSeq: Differential Exon Usage"
    ) +
    base_theme

  # SUPPA2 example
  suppa_data <- data.frame(
    Event = rep(c("Event 1", "Event 2"), each = 2),
    Condition = rep(c("Normoxia", "Hypoxia"), times = 2),
    PSI = c(0.15, 0.8, 0.7, 0.1),
    stringsAsFactors = FALSE
  )
  suppa_data$PSI_percent <- suppa_data$PSI * 100
  suppa_data$Condition <- factor(suppa_data$Condition, levels = c("Normoxia", "Hypoxia"))

  plot_suppa <- ggplot2::ggplot(
    suppa_data,
    ggplot2::aes(x = Event, y = PSI_percent, fill = Condition)
  ) +
    ggplot2::geom_col(position = "dodge", width = 0.75) +
    ggplot2::scale_fill_manual(values = condition_colors) +
    ggplot2::scale_y_continuous(
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100)
    ) +
    ggplot2::labs(
      x = "",
      y = "PSI (%)",
      title = "SUPPA2: Differential Splicing Events"
    ) +
    base_theme

  # DRIMSeq example
  drim_data <- data.frame(
    Transcript = rep(c("Transcript 1", "Transcript 2"), each = 2),
    Condition = rep(c("Normoxia", "Hypoxia"), times = 2),
    Proportion = c(0.7, 0.15, 0.1, 0.8),
    stringsAsFactors = FALSE
  )
  drim_data$Proportion_percent <- drim_data$Proportion * 100
  drim_data$Condition <- factor(drim_data$Condition, levels = c("Normoxia", "Hypoxia"))

  plot_drim <- ggplot2::ggplot(
    drim_data,
    ggplot2::aes(x = Transcript, y = Proportion_percent, fill = Condition)
  ) +
    ggplot2::geom_col(position = "dodge", width = 0.75) +
    ggplot2::scale_fill_manual(values = condition_colors) +
    ggplot2::scale_y_continuous(
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100)
    ) +
    ggplot2::labs(
      x = "",
      y = "Transcript Usage (%)",
      title = "DRIMSeq: Differential Transcript Usage"
    ) +
    base_theme

  grDevices::pdf(
    file.path(outdir, "fig_2_A_splicing_tool_examples.pdf"),
    width = 6,
    height = 5
  )
  print(plot_dex)
  print(plot_suppa)
  print(plot_drim)
  invisible(grDevices::dev.off())

  invisible(list(
    DEXSeq = plot_dex,
    SUPPA2 = plot_suppa,
    DRIMSeq = plot_drim
  ))
}

# --------------------------------------------------
# Figure 2B - Overlap of Splicing Tools by Timepoint
# --------------------------------------------------

plot_splicing_tool_overlap_across_timepoints <- function(
  res,
  outdir = NULL,
  title = "Overlap of Splicing Tools Across Hypoxia Timepoints"
) {
  .check_pkg("eulerr")

  if (is.null(outdir)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_2")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  tool_colors <- unname(splicing_tool_colors[c("DEXSeq", "DRIMSeq", "SUPPA2")])

  dex_by_tp <- res$gene_sets$dex_genes_by_tp
  drim_by_tp <- res$gene_sets$drim_genes_by_tp
  suppa_by_tp <- res$gene_sets$suppa_genes_by_tp

  timepoints <- c("H1", "H3", "H24")
  tp_titles <- c(
    H1 = "Hypoxia (1H)",
    H3 = "Hypoxia (3H)",
    H24 = "Hypoxia (24H)"
  )

  grDevices::pdf(
    file.path(outdir, "fig_2_B_splicing_tool_overlap_across_timepoints.pdf"),
    width = 6,
    height = 6
  )

  for (tp in timepoints) {
    dex_genes <- unique(dex_by_tp[[tp]])
    drim_genes <- unique(drim_by_tp[[tp]])
    suppa_genes <- unique(suppa_by_tp[[tp]])

    fit <- eulerr::euler(
      c(
        "DEXSeq" = length(dex_genes),
        "DRIMSeq" = length(drim_genes),
        "SUPPA2" = length(suppa_genes),
        "DEXSeq&DRIMSeq" = length(intersect(dex_genes, drim_genes)),
        "DEXSeq&SUPPA2" = length(intersect(dex_genes, suppa_genes)),
        "DRIMSeq&SUPPA2" = length(intersect(drim_genes, suppa_genes)),
        "DEXSeq&DRIMSeq&SUPPA2" = length(Reduce(intersect, list(dex_genes, drim_genes, suppa_genes)))
      ),
      input = "union"
    )

    fig2_b <- plot(
      fit,
      fills = tool_colors,
      quantities = TRUE,
      legend = TRUE,
      main = tp_titles[[tp]]
    )

    grid::grid.newpage()
    grid::grid.draw(fig2_b)
  }

  invisible(grDevices::dev.off())
}

# --------------------------------------------------
# Figure 2C - DEG / Splicing Overlap by Tool and Timepoint
# --------------------------------------------------

plot_deg_splicing_overlap_by_tool <- function(
  res,
  base_size = 15,
  outdir = NULL,
  title = "DEG and Splicing Overlap by Tool"
) {
  .check_pkg(c(
    "ggplot2",
    "scales"
  ))

  if (is.null(outdir)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_2")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  df <- as.data.frame(res$overlap$overlap_bar_df, stringsAsFactors = FALSE)

  timepoint_order <- c("H1", "H3", "H24")
  tool_order <- c("DEXSeq", "DRIMSeq", "SUPPA2")

  tp_labels <- c(
    H1  = "Hypoxia (1H)",
    H3  = "Hypoxia (3H)",
    H24 = "Hypoxia (24H)"
  )

  metric_order <- c(
    "% DEGs Also Identified as Alternatively Spliced",
    "% ASGs Also Identified as Differentially Expressed"
  )

  metric_shades <- c(
    "% DEGs Also Identified as Alternatively Spliced" = "DEG_ALSO_AS",
    "% ASGs Also Identified as Differentially Expressed" = "ASG_ALSO_DEG"
  )

  df$timepoint <- factor(df$timepoint, levels = timepoint_order)
  df$method <- factor(df$method, levels = tool_order)
  df$metric <- factor(df$metric, levels = metric_order)

  tool_fill_colors <- .generate_shaded_palette(
    base_colors = splicing_tool_colors[tool_order],
    labels = c("ASG_ALSO_DEG", "BASE", "DEG_ALSO_AS")
  )

  make_one_plot <- function(tp) {
    df_tp <- df[df$timepoint == tp, , drop = FALSE]

    df_tp$fill_key <- paste(df_tp$method, df_tp$metric, sep = "__")

    fill_keys <- as.vector(outer(tool_order, metric_order, paste, sep = "__"))

    fill_values <- setNames(
      vapply(fill_keys, function(key) {
        parts <- strsplit(key, "__", fixed = TRUE)[[1]]
        tool <- parts[1]
        metric <- parts[2]
        shade <- metric_shades[[metric]]

        unname(tool_fill_colors[[tool]][shade])
      }, character(1)),
      fill_keys
    )

    fill_labels <- setNames(
      sub(".*__", "", fill_keys),
      fill_keys
    )

    ggplot2::ggplot(
      df_tp,
      ggplot2::aes(
        x = method,
        y = value,
        fill = fill_key,
        group = metric
      )
    ) +
      ggplot2::geom_col(
        position = ggplot2::position_dodge(width = 0.75),
        width = 0.65
      ) +
      ggplot2::scale_fill_manual(
        values = fill_values,
        breaks = fill_keys,
        labels = fill_labels,
        name = NULL
      ) +
      ggplot2::guides(
        fill = ggplot2::guide_legend(
          nrow = 2,
          byrow = TRUE,
          override.aes = list(
            linewidth = 0,
            color = NA
          )
        )
      ) +
      ggplot2::scale_x_discrete(limits = tool_order) +
      ggplot2::scale_y_continuous(
        labels = scales::percent,
        limits = c(0, 0.52),
        expand = ggplot2::expansion(mult = c(0, 0.03))
      ) +
      ggplot2::labs(
        title = tp_labels[[tp]],
        x = "",
        y = "% Genes"
      ) +
      ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.title = ggplot2::element_blank(),
        legend.key.width = grid::unit(0.35, "cm"),
        legend.key.height = grid::unit(0.35, "cm"),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(vjust = -3)
      )
  }

  fig_h1 <- make_one_plot("H1")
  fig_h3 <- make_one_plot("H3")
  fig_h24 <- make_one_plot("H24")

  grDevices::pdf(
    file.path(outdir, "fig_2_C_deg_splicing_overlap_by_tool.pdf"),
    width = 12,
    height = 9
  )
  print(fig_h1)
  print(fig_h3)
  print(fig_h24)
  invisible(grDevices::dev.off())

  invisible(list(
    H1 = fig_h1,
    H3 = fig_h3,
    H24 = fig_h24
  ))
}

# --------------------------------------------------
# Figure 2D - DEG vs ASG Gene Overlap Across Timepoints
# --------------------------------------------------

plot_splicing_vs_deg_gene_overlap <- function(
  res,
  base_size = 15,
  outdir = NULL,
  title = "Gene-Level Overlap of Differential Expression and Splicing"
) {
  .check_pkg(c(
    "scales",
    "ggplot2",
    "ggnewscale"
  ))

  if (is.null(outdir)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_2")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  tp_labels <- c(
    H1  = "Hypoxia (1H)",
    H3  = "Hypoxia (3H)",
    H24 = "Hypoxia (24H)"
  )

  tp_fill_colors <- .get_timepoint_shaded_colors(
    timepoints = c("H1", "H3", "H24"),
    labels = c("ASG", "OVERLAP", "DEG")
  )

  shade_order <- c("ASG", "OVERLAP", "DEG")
  group_order <- c("ASG Genes", "Shared Genes", "DEG Genes")

  df <- as.data.frame(res$overlap$compare_df, stringsAsFactors = FALSE)
  df$timepoint <- factor(df$timepoint, levels = c("H1", "H3", "H24"))
  df$group <- factor(df$group, levels = group_order)

  df_h1 <- df[df$timepoint == "H1", , drop = FALSE]
  df_h3 <- df[df$timepoint == "H3", , drop = FALSE]
  df_h24 <- df[df$timepoint == "H24", , drop = FALSE]

  df_h1$fill_col <- unname(tp_fill_colors$H1[shade_order][match(df_h1$group, group_order)])
  df_h3$fill_col <- unname(tp_fill_colors$H3[shade_order][match(df_h3$group, group_order)])
  df_h24$fill_col <- unname(tp_fill_colors$H24[shade_order][match(df_h24$group, group_order)])

  fig2_d <- ggplot2::ggplot() +
    ggplot2::geom_col(
      data = df_h1,
      ggplot2::aes(x = timepoint, y = count, fill = fill_col),
      position = ggplot2::position_fill(reverse = TRUE)
    ) +
    ggplot2::scale_fill_identity(
      guide = "legend",
      breaks = unname(tp_fill_colors$H1[shade_order]),
      labels = group_order,
      name = tp_labels["H1"]
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(order = 1)) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_col(
      data = df_h3,
      ggplot2::aes(x = timepoint, y = count, fill = fill_col),
      position = ggplot2::position_fill(reverse = TRUE)
    ) +
    ggplot2::scale_fill_identity(
      guide = "legend",
      breaks = unname(tp_fill_colors$H3[shade_order]),
      labels = group_order,
      name = tp_labels["H3"]
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(order = 2)) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_col(
      data = df_h24,
      ggplot2::aes(x = timepoint, y = count, fill = fill_col),
      position = ggplot2::position_fill(reverse = TRUE)
    ) +
    ggplot2::scale_fill_identity(
      guide = "legend",
      breaks = unname(tp_fill_colors$H24[shade_order]),
      labels = group_order,
      name = tp_labels["H24"]
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(order = 3)) +
    ggplot2::scale_x_discrete(
      limits = c("H1", "H3", "H24"),
      labels = tp_labels
    ) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(
      title = NULL,
      x = "",
      y = "% of Significant Genes"
    ) +
    ggplot2::theme_minimal(base_size = base_size)

  ggplot2::ggsave(
    filename = file.path(outdir, "fig_2_D_splicing_vs_deg_gene_overlap.pdf"),
    plot = fig2_d,
    width = 9,
    height = 5
  )

  fig2_d
}

# --------------------------------------------------
# Convenience Wrapper
# --------------------------------------------------

plot_fig_2_all <- function(
  deseq_results,
  dexseq_results,
  drimseq_results,
  suppa_results,
  outdir = NULL,
  save_tables = FALSE
) {
  if (is.null(outdir)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_2")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  res <- run_transcript_level_analysis(
    deseq_results = deseq_results,
    dexseq_results = dexseq_results,
    drimseq_results = drimseq_results,
    suppa_results = suppa_results,
    outdir = outdir,
    save_tables = save_tables
  )

  p2a <- plot_splicing_tool_examples()
  p2b <- plot_splicing_tool_overlap_across_timepoints(res = res, outdir = outdir)
  p2c <- plot_deg_splicing_overlap_by_tool(res = res, outdir = outdir)
  p2d <- plot_splicing_vs_deg_gene_overlap(res = res, outdir = outdir)

  invisible(list(
    analysis = res,
    fig_2_a = p2a,
    fig_2_b = p2b,
    fig_2_c = p2c,
    fig_2_d = p2d
  ))
}
