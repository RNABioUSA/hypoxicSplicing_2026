# ============================================================
# plot_fig_1.R
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

.scale_rows <- function(x) {
  x <- as.matrix(x)
  m <- apply(x, 1, mean, na.rm = TRUE)
  s <- apply(x, 1, sd, na.rm = TRUE)

  s[is.na(s) | s == 0] <- 1

  out <- (x - m) / s
  out[is.na(out)] <- 0
  out
}

.calc_hm_size <- function(hm, unit = "inch") {
  grDevices::pdf(NULL)
  hm <- ComplexHeatmap::draw(hm)
  w <- ComplexHeatmap:::width(hm)
  w <- grid::convertX(w, unit, valueOnly = TRUE)
  h <- ComplexHeatmap:::height(hm)
  h <- grid::convertY(h, unit, valueOnly = TRUE)
  grDevices::dev.off()
  c(w, h)
}

# --------------------------------------------------
# Gene-Level Analysis
# --------------------------------------------------

run_gene_level_analysis <- function(
  deseq_results,
  n_clusters = 5,
  seed = 333,
  outdir = NULL,
  save_tables = FALSE
) {
  .check_pkg("DESeq2")
  .check_pkg("tibble")
  .check_pkg("dplyr")
  .check_pkg("ComplexHeatmap")

  if (is.null(outdir)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_1")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(deseq_results$dds)) {
    stop("Expected deseq_results$dds.")
  }
  if (is.null(deseq_results$results$gene_sig_by_tp)) {
    stop("Expected deseq_results$results$gene_sig_by_tp.")
  }

  dds <- deseq_results$dds

  sig_genes_by_tp <- .subset_gene_lists_by_tp(
    .get_sig_genes_by_tp(deseq_results$results$gene_sig_by_tp),
    timepoints = c("H1", "H3", "H24")
  )

  sig_deg <- unique(unlist(sig_genes_by_tp, use.names = FALSE))
  sig_deg <- sig_deg[!is.na(sig_deg) & sig_deg != ""]

  counts_normalized <- as.data.frame(DESeq2::counts(dds, normalized = TRUE))
  counts_normalized$ensgene <- .strip_ens_version(rownames(counts_normalized))
  rownames(counts_normalized) <- NULL

  coldata <- as.data.frame(SummarizedExperiment::colData(dds))
  sample_names <- rownames(coldata)

  if (!all(sample_names %in% colnames(counts_normalized))) {
    stop("Sample names in colData(dds) do not match normalized count matrix columns.")
  }
  if (!("condition" %in% colnames(coldata))) {
    stop("Expected 'condition' column in colData(dds).")
  }

  counts_normalized$Normoxia <- rowMeans(
    counts_normalized[, sample_names[coldata$condition == "C1"], drop = FALSE]
  )
  counts_normalized$Hypoxia_1H <- rowMeans(
    counts_normalized[, sample_names[coldata$condition == "H1"], drop = FALSE]
  )
  counts_normalized$Hypoxia_3H <- rowMeans(
    counts_normalized[, sample_names[coldata$condition == "H3"], drop = FALSE]
  )
  counts_normalized$Hypoxia_24H <- rowMeans(
    counts_normalized[, sample_names[coldata$condition == "H24"], drop = FALSE]
  )

  counts_normalized_average <- counts_normalized |>
    dplyr::select(ensgene, Normoxia, Hypoxia_1H, Hypoxia_3H, Hypoxia_24H)

  counts_normalized_average_scaled <- .scale_rows(
    counts_normalized_average[, c("Normoxia", "Hypoxia_1H", "Hypoxia_3H", "Hypoxia_24H")]
  )
  rownames(counts_normalized_average_scaled) <- counts_normalized_average$ensgene

  counts_normalized_average_scaled_sig <- counts_normalized_average_scaled[
    rownames(counts_normalized_average_scaled) %in% sig_deg, ,
    drop = FALSE
  ]

  if (nrow(counts_normalized_average_scaled_sig) == 0) {
    stop("No significant DEGs found for Figure 1.")
  }

  cluster_colors_raw <- c(
    "1" = unname(category_base_colors["Category 1"]),
    "2" = unname(category_base_colors["Category 2"]),
    "3" = unname(category_base_colors["Category 3"]),
    "4" = unname(category_base_colors["Category 4"]),
    "5" = unname(category_base_colors["Category 5"])
  )

  set.seed(seed)
  km <- stats::kmeans(counts_normalized_average_scaled_sig, centers = n_clusters, nstart = 25)

  set.seed(seed)
  tmp_hm <- ComplexHeatmap::Heatmap(
    as.matrix(counts_normalized_average_scaled_sig),
    name = "hm_tmp",
    col = heat_colors,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    cluster_rows = TRUE,
    row_split = km$cluster,
    row_title = NULL,
    show_row_dend = TRUE,
    show_row_names = FALSE,
    show_column_names = FALSE
  )

  grDevices::pdf(NULL)
  tmp_hm <- ComplexHeatmap::draw(tmp_hm)
  cluster_row_order <- ComplexHeatmap::row_order(tmp_hm)
  grDevices::dev.off()

  cluster_display_order <- names(cluster_row_order)

  if (is.null(cluster_display_order) || any(cluster_display_order == "")) {
    cluster_display_order <- as.character(seq_along(cluster_row_order))
  }

  clusters_df <- data.frame(
    ensgene = rownames(counts_normalized_average_scaled_sig),
    cluster = as.character(km$cluster),
    stringsAsFactors = FALSE
  )
  clusters_df <- clusters_df[order(clusters_df$ensgene), , drop = FALSE]
  rownames(clusters_df) <- NULL

  cluster_sizes <- as.data.frame(table(clusters_df$cluster), stringsAsFactors = FALSE)
  colnames(cluster_sizes) <- c("cluster", "n_genes")

  out <- list(
    meta = list(
      outdir = outdir,
      n_clusters = n_clusters,
      seed = seed
    ),
    sig_genes_by_tp = sig_genes_by_tp,
    sig_deg = sig_deg,
    expression = list(
      counts_normalized_average = counts_normalized_average,
      counts_normalized_average_scaled = counts_normalized_average_scaled,
      counts_normalized_average_scaled_sig = counts_normalized_average_scaled_sig
    ),
    clustering = list(
      kmeans = km,
      clusters_df = clusters_df,
      cluster_sizes = cluster_sizes,
      cluster_display_order = cluster_display_order,
      cluster_row_order = cluster_row_order,
      cluster_colors_raw = cluster_colors_raw
    )
  )

  if (isTRUE(save_tables)) {
    utils::write.table(
      clusters_df,
      file = file.path(outdir, "fig_1_clusters.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    utils::write.table(
      as.data.frame(counts_normalized_average_scaled_sig) |>
        tibble::rownames_to_column("ensgene"),
      file = file.path(outdir, "fig_1_scaled_sig_expression.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  saveRDS(out, file.path(outdir, "fig_1_gene_level_analysis.rds"))
  out
}

# --------------------------------------------------
# Figure 1A - DEG Overlap Across Timepoints
# --------------------------------------------------

plot_deg_overlap_across_timepoints <- function(
  res,
  outdir = NULL,
  title = "Overlap of Differentially Expressed Genes Across Hypoxia Timepoints"
) {
  .check_pkg("eulerr")

  if (is.null(outdir)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_1")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  tp_deg_colors <- vapply(
    .get_timepoint_shaded_colors(
      timepoints = c("H1", "H3", "H24"),
      labels = c("TX", "OVERLAP", "GENE")
    ),
    function(x) unname(x["GENE"]),
    character(1)
  )

  genes_by_tp <- res$sig_genes_by_tp
  h1_genes <- unique(genes_by_tp[["H1"]])
  h3_genes <- unique(genes_by_tp[["H3"]])
  h24_genes <- unique(genes_by_tp[["H24"]])

  fit <- eulerr::euler(
    c(
      "1H" = length(h1_genes),
      "3H" = length(h3_genes),
      "24H" = length(h24_genes),
      "1H&3H" = length(intersect(h1_genes, h3_genes)),
      "1H&24H" = length(intersect(h1_genes, h24_genes)),
      "3H&24H" = length(intersect(h3_genes, h24_genes)),
      "1H&3H&24H" = length(Reduce(intersect, list(h1_genes, h3_genes, h24_genes)))
    ),
    input = "union"
  )

  fig1_a <- plot(
    fit,
    fills = unname(tp_deg_colors[c("H1", "H3", "H24")]),
    quantities = TRUE,
    legend = list(
      labels = c("1H", "3H", "24H"),
      title = "Hypoxia"
    ),
    main = NULL
  )

  grid::grid.draw(fig1_a)

  grDevices::pdf(
    file.path(outdir, "fig_1_A_deg_overlap_across_timepoints.pdf"),
    width = 6,
    height = 6
  )
  grid::grid.newpage()
  grid::grid.draw(fig1_a)
  invisible(grDevices::dev.off())

  invisible(fig1_a)
}

# --------------------------------------------------
# Figure 1B - PCA Plot of DEG Clusters
# --------------------------------------------------

plot_deg_cluster_pca <- function(
  res,
  outdir = NULL,
  title = "PCA Plot of Identified Clusters",
  base_size = 15
) {
  .check_pkg("factoextra")
  .check_pkg("ggplot2")

  if (is.null(outdir)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_1")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  km <- res$clustering$kmeans
  expr_mat <- res$expression$counts_normalized_average_scaled_sig

  cluster_colors <- res$clustering$cluster_colors_raw
  cluster_colors_plot <- cluster_colors[sort(unique(as.character(km$cluster)))]

  pca <- stats::prcomp(expr_mat, center = TRUE, scale. = FALSE)
  pve <- (pca$sdev^2) / sum(pca$sdev^2)
  pc1_lab <- sprintf("PC1: %.1f%% Variance", 100 * pve[1])
  pc2_lab <- sprintf("PC2: %.1f%% Variance", 100 * pve[2])

  fig1_b <- suppressMessages(
    suppressWarnings(
      factoextra::fviz_cluster(
        km,
        data = expr_mat,
        stand = FALSE,
        geom = "point",
        ellipse = TRUE,
        ellipse.type = "norm",
        ellipse.level = 0.80,
        shape = 19,
        main = title,
        palette = unname(cluster_colors_plot),
        legend.title = "Cluster",
        ggtheme = ggplot2::theme_classic(base_size = base_size),
        xlab = pc1_lab,
        ylab = pc2_lab
      ) +
        ggplot2::guides(
          fill = ggplot2::guide_legend(
            byrow = TRUE,
            override.aes = list(
              shape = 16,
              fill = NA,
              linetype = 0,
              linewidth = 0,
              stroke = 0,
              alpha = 1,
              size = 3
            )
          )
        ) +
        ggplot2::theme(
          legend.key.height = grid::unit(0.6, "cm"),
          legend.key.width = grid::unit(0.6, "cm"),
          legend.text = ggplot2::element_text(
            margin = ggplot2::margin(t = 4, b = 4, l = 4)
          ),
          legend.title = ggplot2::element_text(
            margin = ggplot2::margin(b = 5)
          )
        )
    )
  )

  suppressMessages(
    suppressWarnings(
      ggplot2::ggsave(
        filename = file.path(outdir, "fig_1_B_deg_cluster_pca.pdf"),
        plot = fig1_b,
        width = 15,
        height = 12,
        pointsize = 0.5
      )
    )
  )

  fig1_b
}

# --------------------------------------------------
# Figure 1C - DEG Cluster Timecourse Plot
# --------------------------------------------------

plot_deg_cluster_timecourse <- function(
  res,
  outdir = NULL,
  title = "Hypoxic Response Over Time",
  base_size = 15
) {
  .check_pkg("ggplot2")
  .check_pkg("dplyr")
  .check_pkg("tibble")
  .check_pkg("tidyr")

  if (is.null(outdir)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_1")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  expr_df <- as.data.frame(res$expression$counts_normalized_average_scaled_sig) |>
    tibble::rownames_to_column("ensgene")

  clusters_df <- res$clustering$clusters_df

  cluster_colors <- res$clustering$cluster_colors_raw
  cluster_colors_plot <- stats::setNames(
    unname(cluster_colors[as.character(1:5)]),
    paste("Cluster", 1:5)
  )

  timecourse <- expr_df |>
    dplyr::left_join(clusters_df, by = "ensgene") |>
    dplyr::mutate(
      cluster_raw = as.character(cluster),
      cluster = factor(
        paste("Cluster", cluster_raw),
        levels = paste("Cluster", 1:5)
      )
    ) |>
    tidyr::pivot_longer(
      cols = c("Normoxia", "Hypoxia_1H", "Hypoxia_3H", "Hypoxia_24H"),
      names_to = "condition",
      values_to = "exp"
    ) |>
    dplyr::mutate(
      time = dplyr::case_when(
        condition == "Normoxia" ~ 0,
        condition == "Hypoxia_1H" ~ 1,
        condition == "Hypoxia_3H" ~ 3,
        condition == "Hypoxia_24H" ~ 24,
        TRUE ~ NA_real_
      )
    )

  fig1_c <- ggplot2::ggplot(
    timecourse,
    ggplot2::aes(x = time, y = exp, color = cluster, fill = cluster)
  ) +
    suppressMessages(
      suppressWarnings(
        ggplot2::geom_smooth(
          method = "loess",
          se = FALSE,
          linewidth = 3,
          linetype = 6
        )
      )
    ) +
    ggplot2::scale_x_continuous(breaks = c(0, 1, 3, 24)) +
    ggplot2::scale_color_manual(values = cluster_colors_plot, labels = 1:5) +
    ggplot2::scale_fill_manual(values = cluster_colors_plot, labels = 1:5) +
    ggplot2::labs(
      title = NULL,
      x = "Time (Hours)",
      y = "Scaled Expression",
      color = "Cluster",
      fill = "Cluster"
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.key.height = grid::unit(0.6, "cm"),
      legend.key.width = grid::unit(0.6, "cm"),
      legend.text = ggplot2::element_text(
        margin = ggplot2::margin(t = 4, b = 4, l = 4)
      ),
      legend.title = ggplot2::element_text(
        margin = ggplot2::margin(b = 5)
      )
    )

  suppressMessages(
    suppressWarnings(
      ggplot2::ggsave(
        filename = file.path(outdir, "fig_1_C_deg_cluster_timecourse.pdf"),
        plot = fig1_c,
        width = 15,
        height = 12
      )
    )
  )

  fig1_c
}

# --------------------------------------------------
# Figure 1D - Heatmap of DEG Clusters
# --------------------------------------------------

plot_deg_cluster_heatmap <- function(
  res,
  outdir = NULL,
  title = "Heatmap of Significant Differentially Expressed Genes"
) {
  .check_pkg("ComplexHeatmap")
  .check_pkg("circlize")

  if (is.null(outdir)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_1")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  expr_mat <- res$expression$counts_normalized_average_scaled_sig
  km <- res$clustering$kmeans

  cluster_colors <- res$clustering$cluster_colors_raw

  cluster_display_order <- res$clustering$cluster_display_order
  row_split <- factor(
    as.character(km$cluster),
    levels = cluster_display_order
  )

  display_timepoint_colors <- timepoint_base_colors[
    c("Normoxia", "Hypoxia (1H)", "Hypoxia (3H)", "Hypoxia (24H)")
  ]

  top_annotation_shades <- .generate_shaded_palette(
    base_colors = display_timepoint_colors,
    labels = c("TX", "OVERLAP", "GENE")
  )

  top_annotation_colors <- vapply(
    top_annotation_shades,
    function(x) unname(x["GENE"]),
    character(1)
  )

  top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    top = ComplexHeatmap::anno_block(
      gp = grid::gpar(fill = unname(top_annotation_colors)),
      labels = names(top_annotation_colors),
      labels_gp = grid::gpar(col = "white", fontsize = 12, fontface = "bold")
    )
  )

  left_annotation <- ComplexHeatmap::rowAnnotation(
    left = ComplexHeatmap::anno_block(
      gp = grid::gpar(fill = unname(cluster_colors[cluster_display_order])),
      labels = paste("Cluster", cluster_display_order),
      labels_gp = grid::gpar(col = "white", fontsize = 12, fontface = "bold")
    )
  )

  set.seed(res$meta$seed)
  hm <- ComplexHeatmap::Heatmap(
    as.matrix(expr_mat),
    name = "hm",
    col = rev(heat_colors),
    top_annotation = top_annotation,
    left_annotation = left_annotation,
    heatmap_legend_param = list(
      title = " ",
      legend_height = grid::unit(6, "cm")
    ),
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    cluster_rows = TRUE,
    cluster_row_slices = FALSE,
    column_split = 1:4,
    column_gap = grid::unit(0, "mm"),
    column_title = NULL,
    row_dend_width = grid::unit(1.5, "cm"),
    row_gap = grid::unit(3, "mm"),
    row_split = row_split,
    row_title = rep(" ", 5),
    show_column_names = FALSE,
    show_row_dend = TRUE,
    show_row_names = FALSE,
    border = FALSE,
    show_parent_dend_line = FALSE,
    width = ncol(expr_mat) * grid::unit(1.5, "in"),
    height = grid::unit(12, "in")
  )

  hm_size <- .calc_hm_size(hm)

  grDevices::pdf(
    file.path(outdir, "fig_1_D_deg_cluster_heatmap.pdf"),
    width = hm_size[1],
    height = hm_size[2]
  )

  hm <- ComplexHeatmap::draw(hm)

  for (i in 1:5) {
    for (j in 1:4) {
      ComplexHeatmap::decorate_heatmap_body("hm",
        {
          grid::grid.lines(
            c(0, 1),
            grid::unit(c(0, 0), "native"),
            gp = grid::gpar(lty = 1, lwd = 2.5, col = "black")
          )
        },
        row_slice = i,
        column_slice = j
      )

      ComplexHeatmap::decorate_heatmap_body("hm",
        {
          grid::grid.lines(
            c(0, 1),
            grid::unit(c(1, 1), "native"),
            gp = grid::gpar(lty = 1, lwd = 2.5, col = "black")
          )
        },
        row_slice = i,
        column_slice = j
      )
    }
  }

  for (i in 1:5) {
    ComplexHeatmap::decorate_heatmap_body("hm",
      {
        grid::grid.lines(
          c(0, 0),
          gp = grid::gpar(lty = 1, lwd = 1.75, col = "black")
        )
      },
      row_slice = i,
      column_slice = 1
    )

    ComplexHeatmap::decorate_heatmap_body("hm",
      {
        grid::grid.lines(
          c(1, 1),
          gp = grid::gpar(lty = 1, lwd = 2.5, col = "black")
        )
      },
      row_slice = i,
      column_slice = 4
    )
  }

  invisible(grDevices::dev.off())
  hm
}

# --------------------------------------------------
# Convenience Wrapper
# --------------------------------------------------

plot_fig_1_all <- function(
  deseq_results,
  n_clusters = 5,
  seed = 333,
  outdir = NULL,
  save_tables = FALSE
) {
  if (is.null(outdir)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_1")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  res <- run_gene_level_analysis(
    deseq_results = deseq_results,
    n_clusters = n_clusters,
    seed = seed,
    outdir = outdir,
    save_tables = save_tables
  )

  fig1_a <- plot_deg_overlap_across_timepoints(res = res, outdir = outdir)
  fig1_b <- plot_deg_cluster_pca(res = res, outdir = outdir)
  fig1_c <- plot_deg_cluster_timecourse(res = res, outdir = outdir)
  fig1_d <- plot_deg_cluster_heatmap(res = res, outdir = outdir)

  invisible(list(
    analysis = res,
    fig1_a = fig1_a,
    fig1_b = fig1_b,
    fig1_c = fig1_c,
    fig1_d = fig1_d
  ))
}
