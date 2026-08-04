# ============================================================
# plot_fig_5.R
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
# Transcript Level Enrichment Analysis
# -------------------------

run_splicing_pathway_enrichment_analysis <- function(dexseq_results,
                                                     drimseq_results,
                                                     suppa_results,
                                                     deseq_results,
                                                     timepoints = c("H1", "H3", "H24"),
                                                     enrichment_cutoff = 0.10,
                                                     annot_dataset = "GO:0008150",
                                                     outdir = NULL,
                                                     save_tables = FALSE,
                                                     force_enrich = FALSE) {
  if (is.null(outdir)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_5")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  annot_df <- .get_annot()

  tested_universe_by_tp <- .get_active_gene_universe_by_tp(
    deseq_results = deseq_results,
    timepoints = timepoints
  )

  dex_genes_by_tp <- .get_sig_genes_by_tp(dexseq_results$results$gene_sig_by_tp)
  drim_genes_by_tp <- .get_sig_genes_by_tp(drimseq_results$results$gene_sig_by_tp)
  suppa_genes_by_tp <- .get_sig_genes_by_tp(suppa_results$results$gene_sig_by_tp)
  deg_genes_by_tp <- .get_sig_genes_by_tp(deseq_results$results$gene_sig_by_tp)

  all_splicing_by_tp <- .union_gene_lists_by_tp(
    dex_genes_by_tp,
    drim_genes_by_tp,
    suppa_genes_by_tp
  )

  as_pathways_gene_set_by_tp <- .subtract_gene_lists_by_tp(all_splicing_by_tp, deg_genes_by_tp)
  shared_gene_set_by_tp <- .intersect_gene_lists_by_tp(all_splicing_by_tp, deg_genes_by_tp)
  deg_only_by_tp <- .subtract_gene_lists_by_tp(deg_genes_by_tp, all_splicing_by_tp)

  all_splicing_by_tp <- .filter_gene_lists_to_universe_by_tp(all_splicing_by_tp, tested_universe_by_tp)
  as_pathways_gene_set_by_tp <- .filter_gene_lists_to_universe_by_tp(as_pathways_gene_set_by_tp, tested_universe_by_tp)
  shared_gene_set_by_tp <- .filter_gene_lists_to_universe_by_tp(shared_gene_set_by_tp, tested_universe_by_tp)
  deg_only_by_tp <- .filter_gene_lists_to_universe_by_tp(deg_only_by_tp, tested_universe_by_tp)
  deg_genes_by_tp <- .filter_gene_lists_to_universe_by_tp(deg_genes_by_tp, tested_universe_by_tp)

  cache_dir <- file.path(outdir, "cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  as_pathways_enrich_by_tp <- stats::setNames(lapply(timepoints, function(tp) {
    cache_file <- file.path(
      cache_dir,
      paste0("as_pathways_", tp, "_", gsub("[: ]", "_", annot_dataset), "_fdr", enrichment_cutoff, ".rds")
    )

    .run_panther_enrich_from_ensgenes(
      query_ensgenes = as_pathways_gene_set_by_tp[[tp]],
      universe_ensgenes = tested_universe_by_tp[[tp]],
      annot_df = annot_df,
      organism = 9606,
      ref_organism = 9606,
      annot_dataset = annot_dataset,
      cutoff = enrichment_cutoff,
      cache_file = cache_file,
      force = force_enrich
    )
  }), timepoints)

  deg_pathways_enrich_by_tp <- stats::setNames(lapply(timepoints, function(tp) {
    cache_file <- file.path(
      cache_dir,
      paste0("deg_pathways_", tp, "_", gsub("[: ]", "_", annot_dataset), "_fdr", enrichment_cutoff, ".rds")
    )

    .run_panther_enrich_from_ensgenes(
      query_ensgenes = deg_genes_by_tp[[tp]],
      universe_ensgenes = tested_universe_by_tp[[tp]],
      annot_df = annot_df,
      organism = 9606,
      ref_organism = 9606,
      annot_dataset = annot_dataset,
      cutoff = enrichment_cutoff,
      cache_file = cache_file,
      force = force_enrich
    )
  }), timepoints)

  as_pathways_list <- Filter(Negate(is.null), lapply(timepoints, function(tp) {
    x <- as_pathways_enrich_by_tp[[tp]]$result
    if (is.null(x) || nrow(x) == 0) {
      return(NULL)
    }
    x$timepoint <- tp
    x
  }))

  deg_pathways_list <- Filter(Negate(is.null), lapply(timepoints, function(tp) {
    x <- deg_pathways_enrich_by_tp[[tp]]$result
    if (is.null(x) || nrow(x) == 0) {
      return(NULL)
    }
    x$timepoint <- tp
    x
  }))

  as_pathways_tbl <- if (length(as_pathways_list) == 0) data.frame() else do.call(rbind, as_pathways_list)
  deg_pathways_tbl <- if (length(deg_pathways_list) == 0) data.frame() else do.call(rbind, deg_pathways_list)

  as_pathways_terms_by_tp <- lapply(as_pathways_enrich_by_tp, function(x) {
    res_tp <- x$result
    if (is.null(res_tp) || nrow(res_tp) == 0) character(0) else unique(as.character(res_tp$term.id))
  })

  h1_terms <- unique(as_pathways_terms_by_tp[["H1"]])
  h3_terms <- unique(as_pathways_terms_by_tp[["H3"]])
  h24_terms <- unique(as_pathways_terms_by_tp[["H24"]])

  as_pathways_overlap <- list(
    h1_only = setdiff(h1_terms, union(h3_terms, h24_terms)),
    h3_only = setdiff(h3_terms, union(h1_terms, h24_terms)),
    h24_only = setdiff(h24_terms, union(h1_terms, h3_terms)),
    h1_h3 = setdiff(intersect(h1_terms, h3_terms), h24_terms),
    h1_h24 = setdiff(intersect(h1_terms, h24_terms), h3_terms),
    h3_h24 = setdiff(intersect(h3_terms, h24_terms), h1_terms),
    h1_h3_h24 = Reduce(intersect, list(h1_terms, h3_terms, h24_terms))
  )

  as_pathways_overlap_df <- data.frame(
    group = c("H1 only", "H3 only", "H24 only", "H1 & H3", "H1 & H24", "H3 & H24", "H1 & H3 & H24"),
    count = c(
      length(as_pathways_overlap$h1_only),
      length(as_pathways_overlap$h3_only),
      length(as_pathways_overlap$h24_only),
      length(as_pathways_overlap$h1_h3),
      length(as_pathways_overlap$h1_h24),
      length(as_pathways_overlap$h3_h24),
      length(as_pathways_overlap$h1_h3_h24)
    ),
    stringsAsFactors = FALSE
  )

  pathway_compare_df <- do.call(rbind, lapply(timepoints, function(tp) {
    as_res <- as_pathways_enrich_by_tp[[tp]]$result
    deg_res <- deg_pathways_enrich_by_tp[[tp]]$result

    as_terms <- if (is.null(as_res) || nrow(as_res) == 0) {
      character(0)
    } else {
      unique(as.character(as_res$term.id))
    }

    deg_terms <- if (is.null(deg_res) || nrow(deg_res) == 0) {
      character(0)
    } else {
      unique(as.character(deg_res$term.id))
    }

    data.frame(
      timepoint = tp,
      group = c("ASG Pathways", "Shared Pathways", "DEG Pathways"),
      count = c(
        length(setdiff(as_terms, deg_terms)),
        length(intersect(as_terms, deg_terms)),
        length(setdiff(deg_terms, as_terms))
      ),
      stringsAsFactors = FALSE
    )
  }))

  out <- list(
    meta = list(
      timepoints = timepoints,
      enrichment_cutoff = enrichment_cutoff,
      annot_dataset = annot_dataset,
      outdir = outdir
    ),
    gene_sets = list(
      dex_genes_by_tp = dex_genes_by_tp,
      drim_genes_by_tp = drim_genes_by_tp,
      suppa_genes_by_tp = suppa_genes_by_tp,
      deg_genes_by_tp = deg_genes_by_tp,
      all_splicing_by_tp = all_splicing_by_tp,
      as_pathways_gene_set_by_tp = as_pathways_gene_set_by_tp,
      shared_gene_set_by_tp = shared_gene_set_by_tp,
      deg_only_by_tp = deg_only_by_tp,
      tested_universe_by_tp = tested_universe_by_tp
    ),
    annotation = list(annot_df = annot_df),
    enrichment = list(
      as_pathways_enrich_by_tp = as_pathways_enrich_by_tp,
      deg_pathways_enrich_by_tp = deg_pathways_enrich_by_tp,
      as_pathways_tbl = as_pathways_tbl,
      deg_pathways_tbl = deg_pathways_tbl
    ),
    overlap = list(
      as_pathways_terms_by_tp = as_pathways_terms_by_tp,
      as_pathways_overlap = as_pathways_overlap,
      as_pathways_overlap_df = as_pathways_overlap_df,
      pathway_compare_df = pathway_compare_df
    )
  )

  if (isTRUE(save_tables)) {
    if (nrow(as_pathways_tbl) > 0) utils::write.table(as_pathways_tbl, file.path(outdir, "as_pathways.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    if (nrow(deg_pathways_tbl) > 0) utils::write.table(deg_pathways_tbl, file.path(outdir, "deg_pathways.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(as_pathways_overlap_df, file.path(outdir, "as_pathways_overlap_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(pathway_compare_df, file.path(outdir, "as_vs_deg_pathway_overlap_by_timepoint.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  saveRDS(out, file.path(outdir, "fig_5_splicing_pathway_enrichment_analysis.rds"))
  out
}

# -------------------------
# Pathway Representation Analysis
# -------------------------

run_splicing_pathway_representation_analysis <- function(dexseq_results,
                                                         drimseq_results,
                                                         suppa_results,
                                                         deseq_results,
                                                         annot_datasets = c("GO:0008150", "GO:0003674"),
                                                         pathway_info = NULL,
                                                         outdir = NULL,
                                                         resources_dir = NULL,
                                                         save_tables = FALSE,
                                                         force_mapping = FALSE,
                                                         force_obo = FALSE,
                                                         force_term2gene = FALSE) {
  if (is.null(outdir)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_5")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(resources_dir)) {
    resources_dir <- file.path(PROJECT_ROOT, "resources")
  }

  if (is.null(pathway_info)) {
    pathway_info <- data.frame(
      pathway_id = c(
        "GO:0071456", "GO:0032452", "GO:0006325",
        "GO:0000278", "GO:0006260", "GO:0006284",
        "GO:0001525", "GO:0001936", "GO:0085029",
        "GO:0006412", "GO:0000398", "GO:0033108",
        "GO:0030097", "GO:0090279", "GO:0010594"
      ),
      pathway_label = c(
        "Hypoxic Response",
        "Histone Demethylation",
        "Chromatin Organization",
        "Cell Cycle",
        "DNA Replication",
        "DNA Repair (BER)",
        "Angiogenesis",
        "Endothelial Cell Proliferation",
        "Extracellular Matrix Assembly",
        "Translation",
        "Alternative Splicing",
        "Mitochondrial Bioenergetics",
        "Hematopoiesis",
        "Calcium Ion Import",
        "Cell Migration / Motility"
      ),
      cluster = c(
        "Cluster 1", "Cluster 1", "Cluster 1",
        "Cluster 2", "Cluster 2", "Cluster 2",
        "Cluster 3", "Cluster 3", "Cluster 3",
        "Cluster 4", "Cluster 4", "Cluster 4",
        "Cluster 5", "Cluster 5", "Cluster 5"
      ),
      stringsAsFactors = FALSE
    )
  }
  pathway_info$pathway_label <- as.character(pathway_info$pathway_label)

  annot_df <- .get_annot()

  tested_universe_by_tp <- .get_active_gene_universe_by_tp(
    deseq_results = deseq_results,
    timepoints = c("H1", "H3", "H24")
  )

  dex_genes_by_tp <- .get_sig_genes_by_tp(dexseq_results$results$gene_sig_by_tp)
  drim_genes_by_tp <- .get_sig_genes_by_tp(drimseq_results$results$gene_sig_by_tp)
  suppa_genes_by_tp <- .get_sig_genes_by_tp(suppa_results$results$gene_sig_by_tp)
  deg_genes_by_tp <- .get_sig_genes_by_tp(deseq_results$results$gene_sig_by_tp)

  all_splicing_by_tp <- .union_gene_lists_by_tp(dex_genes_by_tp, drim_genes_by_tp, suppa_genes_by_tp)
  as_pathways_gene_set_by_tp <- .subtract_gene_lists_by_tp(all_splicing_by_tp, deg_genes_by_tp)

  all_splicing_by_tp <- .filter_gene_lists_to_universe_by_tp(all_splicing_by_tp, tested_universe_by_tp)
  as_pathways_gene_set_by_tp <- .filter_gene_lists_to_universe_by_tp(as_pathways_gene_set_by_tp, tested_universe_by_tp)
  deg_genes_by_tp <- .filter_gene_lists_to_universe_by_tp(deg_genes_by_tp, tested_universe_by_tp)

  res_stub <- list(
    meta = list(outdir = outdir),
    gene_sets = list(
      tested_universe_by_tp = tested_universe_by_tp,
      deg_genes_by_tp = deg_genes_by_tp,
      as_pathways_gene_set_by_tp = as_pathways_gene_set_by_tp
    ),
    annotation = list(
      annot_df = annot_df
    )
  )

  mapping_res <- .get_panther_mapping(
    res = res_stub,
    annot_datasets = annot_datasets,
    annot_df = annot_df,
    resources_dir = resources_dir,
    force = force_mapping,
    save_tables = save_tables
  )

  term2gene_res <- .get_term2gene(
    mapping_res = mapping_res,
    resources_dir = resources_dir,
    force_obo = force_obo,
    force = force_term2gene,
    save_tables = save_tables
  )

  term2gene_long <- term2gene_res$long
  term2gene_long <- term2gene_long[
    term2gene_long$pathway_id %in% pathway_info$pathway_id,
    c("pathway_id", "pathway_label", "symbol"),
    drop = FALSE
  ]

  # force display labels from pathway_info
  label_map <- stats::setNames(pathway_info$pathway_label, pathway_info$pathway_id)
  term2gene_long$pathway_label <- unname(label_map[term2gene_long$pathway_id])

  cluster_pathways <- split(term2gene_long["symbol"], term2gene_long$pathway_id)
  cluster_pathways <- lapply(cluster_pathways, function(df) {
    df <- unique(df)
    df <- df[!is.na(df$symbol) & df$symbol != "", , drop = FALSE]
    rownames(df) <- NULL
    df
  })

  deg_ens <- unique(unlist(deg_genes_by_tp, use.names = FALSE))
  as_ens <- unique(unlist(as_pathways_gene_set_by_tp, use.names = FALSE))

  deg_symbols <- annot_df$symbol[match(deg_ens, annot_df$ensgene)]
  deg_symbols <- unique(stats::na.omit(as.character(deg_symbols)))

  as_symbols <- annot_df$symbol[match(as_ens, annot_df$ensgene)]
  as_symbols <- unique(stats::na.omit(as.character(as_symbols)))

  wide_list <- lapply(names(cluster_pathways), function(pid) {
    pathway_df <- cluster_pathways[[pid]]

    info_row <- pathway_info[pathway_info$pathway_id == pid, , drop = FALSE]
    if (nrow(info_row) != 1) {
      stop("Expected exactly one row in pathway_info for pathway_id: ", pid)
    }

    pathway_symbols <- unique(stats::na.omit(as.character(pathway_df$symbol)))
    pathway_size <- length(pathway_symbols)

    deg_hits <- intersect(pathway_symbols, deg_symbols)
    as_hits <- intersect(pathway_symbols, as_symbols)

    shared_hits <- intersect(deg_hits, as_hits)
    as_only_hits <- setdiff(as_hits, deg_hits)
    deg_only_hits <- setdiff(deg_hits, as_hits)

    data.frame(
      pathway_id = pid,
      pathway_label = info_row$pathway_label,
      cluster = info_row$cluster,
      pathway_size = pathway_size,
      n_as_only = length(as_only_hits),
      n_shared = length(shared_hits),
      n_deg_only = length(deg_only_hits),
      p_as_only = if (pathway_size > 0) length(as_only_hits) / pathway_size else NA_real_,
      p_shared = if (pathway_size > 0) length(shared_hits) / pathway_size else NA_real_,
      p_deg_only = if (pathway_size > 0) length(deg_only_hits) / pathway_size else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  wide_df <- do.call(rbind, wide_list)

  .make_long <- function(wide_df, value_cols, group_labels) {
    do.call(rbind, lapply(seq_along(value_cols), function(i) {
      data.frame(
        pathway_id = wide_df$pathway_id,
        pathway_label = wide_df$pathway_label,
        cluster = wide_df$cluster,
        pathway_size = wide_df$pathway_size,
        group = group_labels[i],
        value = wide_df[[value_cols[i]]],
        stringsAsFactors = FALSE
      )
    }))
  }

  long_perc <- .make_long(
    wide_df,
    value_cols = c("p_deg_only", "p_shared", "p_as_only"),
    group_labels = c("DEG Pathways", "Shared Pathways", "ASG Pathways")
  )

  pathway_representation <- list(
    wide = wide_df,
    long_perc = long_perc
  )

  out <- list(
    meta = list(outdir = outdir),
    pathway_info = pathway_info,
    mapping = mapping_res,
    term2gene = term2gene_res,
    cluster_pathways = cluster_pathways,
    pathway_representation = pathway_representation
  )

  if (isTRUE(save_tables)) {
    utils::write.table(
      pathway_representation$wide,
      file = file.path(outdir, "clustered_pathway_representation.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  saveRDS(out, file.path(outdir, "fig_5_splicing_pathway_representation_analysis.rds"))
  out
}

# --------------------------------------------------
# Figure 5A - Plot Pathway Representation
# --------------------------------------------------

plot_clustered_pathway_representation <- function(
  representation_res,
  pathway_info = NULL,
  cluster_palette = NULL,
  base_size = 15,
  outdir = NULL,
  title = "Gene- and Transcript-Level Representation of Selected Hypoxia Pathways"
) {
  .check_pkg("ggplot2")
  .check_pkg("scales")
  .check_pkg("ggnewscale")

  if (is.null(outdir)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_5")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(pathway_info)) {
    pathway_info <- representation_res$pathway_info
  }

  pathway_representation <- representation_res$pathway_representation

  cluster_order <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5")

  # legend order
  shade_order <- c("ASG", "OVERLAP", "DEG")

  # enforce cluster order
  pathway_info$cluster <- factor(pathway_info$cluster, levels = cluster_order)

  # preserve row order exactly as defined
  pathway_info$.row_id <- seq_len(nrow(pathway_info))
  pathway_info <- pathway_info[order(pathway_info$cluster, pathway_info$.row_id), ]

  # keep the order exactly as defined
  pathway_levels <- pathway_info$pathway_label

  # cluster colors
  if (is.null(cluster_palette)) {
    cluster_base_colors <- c(
      "Cluster 1" = unname(category_base_colors["Category 1"]),
      "Cluster 2" = unname(category_base_colors["Category 2"]),
      "Cluster 3" = unname(category_base_colors["Category 3"]),
      "Cluster 4" = unname(category_base_colors["Category 4"]),
      "Cluster 5" = unname(category_base_colors["Category 5"])
    )
  } else {
    cluster_base_colors <- cluster_palette
  }

  cluster_base_colors <- cluster_base_colors[cluster_order]

  cluster_signal_colors <- .generate_shaded_palette(
    cluster_base_colors,
    labels = c("ASG", "OVERLAP", "DEG")
  )

  df <- pathway_representation$long_perc
  wide_df <- pathway_representation$wide

  df$pathway_label <- factor(df$pathway_label, levels = pathway_levels)
  wide_df$pathway_label <- factor(wide_df$pathway_label, levels = pathway_levels)

  df$cluster <- factor(df$cluster, levels = cluster_order)
  wide_df$cluster <- factor(wide_df$cluster, levels = cluster_order)

  # map groups
  df$fill_key <- NA
  df$fill_key[df$group == "ASG Pathways"] <- "ASG"
  df$fill_key[df$group == "Shared Pathways"] <- "OVERLAP"
  df$fill_key[df$group == "DEG Pathways"] <- "DEG"

  df$fill_key <- factor(df$fill_key, levels = shade_order)

  # annotation rows
  row_df <- unique(wide_df[, c("pathway_label", "cluster")])

  # use the displayed y-axis order
  display_levels <- rev(pathway_levels)
  row_df$row_id <- match(as.character(row_df$pathway_label), display_levels)

  ann_df <- aggregate(
    row_id ~ cluster,
    data = row_df,
    FUN = function(x) {
      c(
        ymin = min(x) - 0.35,
        ymax = max(x) + 0.35,
        y = mean(range(x))
      )
    }
  )

  ann_df <- do.call(data.frame, ann_df)
  colnames(ann_df) <- c("cluster", "ymin", "ymax", "y")

  ann_df$cluster <- factor(as.character(ann_df$cluster), levels = cluster_order)
  ann_df <- ann_df[order(ann_df$cluster), ]

  xmin_lim <- -0.09
  block_width <- 0.075

  ann_df$xmin <- xmin_lim
  ann_df$xmax <- xmin_lim + block_width
  ann_df$fill <- unname(cluster_base_colors[as.character(ann_df$cluster)])

  fig5_a <- ggplot2::ggplot()

  for (i in seq_along(cluster_order)) {
    cl <- cluster_order[i]
    df_cl <- df[df$cluster == cl, ]

    fig5_a <- fig5_a +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_col(
        data = df_cl,
        ggplot2::aes(x = value, y = pathway_label, fill = fill_key),
        width = 0.75,
        position = "stack"
      ) +
      ggplot2::scale_fill_manual(
        values = unname(cluster_signal_colors[[cl]][shade_order]),
        breaks = shade_order,
        labels = c("ASG", "OVERLAP", "DEG"),
        name = cl,
        guide = ggplot2::guide_legend(order = i)
      )
  }

  fig5_a <- fig5_a +
    ggplot2::geom_rect(
      data = ann_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = ann_df$fill,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = ann_df,
      ggplot2::aes(x = (xmin + xmax) / 2, y = y, label = cluster),
      angle = 90,
      fontface = "bold",
      color = "white",
      inherit.aes = FALSE
    ) +
    ggplot2::scale_x_continuous(
      limits = c(xmin_lim, 1),
      labels = scales::percent
    ) +
    ggplot2::scale_y_discrete(limits = rev(pathway_levels)) +
    ggplot2::labs(
      title = NULL,
      x = "Pathway Gene Representation (%)",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(hjust = 0),
      panel.spacing.y = grid::unit(0.6, "lines"),
      axis.text.y = ggplot2::element_text(size = base_size - 1)
    )

  ggplot2::ggsave(
    filename = file.path(outdir, "fig_5_A_clustered_pathway_representation.pdf"),
    plot = fig5_a,
    width = 10,
    height = 10
  )

  return(fig5_a)
}

# --------------------------------------------------
# Figure 5B -  Plot AS Pathway Overlap Across Timepoints
# --------------------------------------------------

plot_splicing_pathway_overlap_across_timepoints <- function(
  res,
  outdir = NULL,
  title = "Transcript-Level (ASG) Pathways"
) {
  .check_pkg("eulerr")

  if (is.null(outdir)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_5")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  tp_as_colors <- vapply(
    .get_timepoint_shaded_colors(
      timepoints = c("H1", "H3", "H24"),
      labels = c("ASG", "OVERLAP", "DEG")
    ),
    function(x) unname(x["ASG"]),
    character(1)
  )

  terms_by_tp <- res$overlap$as_pathways_terms_by_tp
  h1_terms <- unique(terms_by_tp[["H1"]])
  h3_terms <- unique(terms_by_tp[["H3"]])
  h24_terms <- unique(terms_by_tp[["H24"]])

  fit <- eulerr::euler(
    c(
      "1H" = length(h1_terms),
      "3H" = length(h3_terms),
      "24H" = length(h24_terms),
      "1H&3H" = length(intersect(h1_terms, h3_terms)),
      "1H&24H" = length(intersect(h1_terms, h24_terms)),
      "3H&24H" = length(intersect(h3_terms, h24_terms)),
      "1H&3H&24H" = length(Reduce(intersect, list(h1_terms, h3_terms, h24_terms)))
    ),
    input = "union"
  )

  fig5_b <- plot(
    fit,
    fills = unname(tp_as_colors[c("H1", "H3", "H24")]),
    quantities = TRUE,
    legend = list(
      labels = c("1H", "3H", "24H"),
      title = "Hypoxia"
    ),
    main = title
  )

  grid::grid.draw(fig5_b)

  pdf(file.path(outdir, "fig_5_B_splicing_pathway_overlap_across_timepoints.pdf"), width = 6, height = 6)
  grid::grid.newpage()
  grid::grid.draw(fig5_b)
  invisible(dev.off())
}

# --------------------------------------------------
# Figure 5C - Plot AS vs DEG Pathway Overlap
# --------------------------------------------------

plot_splicing_vs_deg_pathway_overlap <- function(
  res,
  base_size = 15,
  outdir = NULL,
  title = "Gene- and Transcript-Level Pathway Enrichment"
) {
  .check_pkg("ggplot2")
  .check_pkg("scales")
  .check_pkg("ggnewscale")

  if (is.null(outdir)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_5")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  tp_labels <- c(
    H1 = "1H",
    H3 = "3H",
    H24 = "24H"
  )

  tp_fill_colors <- .get_timepoint_shaded_colors(
    timepoints = c("H1", "H3", "H24"),
    labels = c("ASG", "OVERLAP", "DEG")
  )

  shade_order <- c("ASG", "OVERLAP", "DEG")
  group_order <- c("ASG Pathways", "Shared Pathways", "DEG Pathways")

  df <- as.data.frame(res$overlap$pathway_compare_df, stringsAsFactors = FALSE)
  df$timepoint <- factor(df$timepoint, levels = c("H1", "H3", "H24"))
  df$group <- factor(df$group, levels = group_order)

  df_h1 <- df[df$timepoint == "H1", , drop = FALSE]
  df_h3 <- df[df$timepoint == "H3", , drop = FALSE]
  df_h24 <- df[df$timepoint == "H24", , drop = FALSE]

  df_h1$fill_col <- unname(tp_fill_colors$H1[shade_order][match(df_h1$group, group_order)])
  df_h3$fill_col <- unname(tp_fill_colors$H3[shade_order][match(df_h3$group, group_order)])
  df_h24$fill_col <- unname(tp_fill_colors$H24[shade_order][match(df_h24$group, group_order)])

  fig5_c <- ggplot2::ggplot() +
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
      labels = c(
        H1  = "1H",
        H3  = "3H",
        H24 = "24H"
      )
    ) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(
      title = NULL,
      x = "Hypoxia",
      y = "% of Enriched Pathways"
    ) +
    ggplot2::theme_minimal(base_size = base_size)

  ggplot2::ggsave(
    filename = file.path(outdir, "fig_5_C_splicing_vs_deg_pathway_overlap.pdf"),
    plot = fig5_c,
    width = 9,
    height = 5
  )

  return(fig5_c)
}

# --------------------------------------------------
# Figure 5D - Plot Unique AS Pathways
# --------------------------------------------------

.selected_as_pathway_terms <- function() {
  list(
    H1 = c(
      "GO:0001568", # blood vessel development
      "GO:0001525", # angiogenesis
      "GO:0043542", # endothelial cell migration
      "GO:2000145", # regulation of cell motility
      "GO:0048646", # anatomical structure formation involved in morphogenesis
      "GO:0043067", # regulation of programmed cell death
      "GO:0042981", # regulation of apoptotic process
      "GO:0110076", # negative regulation of ferroptosis
      "GO:1902806", # regulation of cell cycle G1/S phase transition
      "GO:0072331", # signal transduction by p53 class mediator
      "GO:0002181", # cytoplasmic translation
      "GO:0006325" # chromatin organization
    ),
    H3 = c(
      "GO:0008380", # RNA splicing
      "GO:0006397", # mRNA processing
      "GO:0006914", # autophagy
      "GO:0016236", # macroautophagy
      "GO:0000045", # autophagosome assembly
      "GO:0051169", # nuclear transport
      "GO:0048193", # Golgi vesicle transport
      "GO:0016192", # vesicle-mediated transport
      "GO:0006405", # RNA export from nucleus
      "GO:0006974", # DNA damage response
      "GO:0033365", # protein localization to organelle
      "GO:0045184" # establishment of protein localization
    ),
    H24 = c(
      "GO:1900037", # regulation of cellular response to hypoxia
      "GO:1900038", # negative regulation of cellular response to hypoxia
      "GO:0060765", # regulation of androgen receptor signaling pathway
      "GO:2001028", # positive regulation of endothelial cell chemotaxis
      "GO:0034063", # stress granule assembly
      "GO:0006397", # mRNA processing
      "GO:0006413", # translational initiation
      "GO:0006446", # regulation of translational initiation
      "GO:0002183", # cytoplasmic translational initiation
      "GO:0072599", # establishment of protein localization to endoplasmic reticulum
      "GO:0045047", # protein targeting to ER
      "GO:0006892" # post-Golgi vesicle-mediated transport
    )
  )
}

plot_transcript_unique_pathway_dotplot_curated <- function(
  res,
  selected_terms = .selected_as_pathway_terms(),
  base_size = 15,
  wrap_width = 28,
  outdir = NULL,
  filename = "fig_5_D_transcript_unique_pathways_curated_dotplot.pdf"
) {
  .check_pkg(c("dplyr", "ggplot2", "stringr", "scales", "ggh4x"))

  if (is.null(outdir)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_5")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  dot_df <- dplyr::bind_rows(lapply(names(selected_terms), function(tp) {
    as_res <- res$enrichment$as_pathways_enrich_by_tp[[tp]]$result
    deg_res <- res$enrichment$deg_pathways_enrich_by_tp[[tp]]$result
    if (is.null(as_res) || nrow(as_res) == 0) {
      return(NULL)
    }

    deg_terms <- if (is.null(deg_res) || nrow(deg_res) == 0) character(0) else unique(as.character(deg_res$term.id))
    keep <- selected_terms[[tp]]

    as_res |>
      dplyr::filter(term.id %in% keep, !(term.id %in% deg_terms)) |>
      dplyr::mutate(
        timepoint = tp,
        term_order = match(term.id, keep),
        gene_ratio = number_in_list / number_in_reference,
        term_clean = term.label |>
          gsub("_", " ", x = _) |>
          stringr::str_to_lower() |>
          stringr::str_to_title() |>
          gsub("\\bRna\\b", "RNA", x = _) |>
          gsub("\\bDna\\b", "DNA", x = _) |>
          gsub("\\bMrna\\b", "mRNA", x = _) |>
          gsub("\\bEr\\b", "ER", x = _) |>
          gsub("\\bAtp\\b", "ATP", x = _) |>
          gsub("\\bAdp\\b", "ADP", x = _) |>
          gsub("\\bAmp\\b", "AMP", x = _) |>
          stringr::str_wrap(width = wrap_width)
      )
  }))

  if (is.null(dot_df) || nrow(dot_df) == 0) {
    stop("No selected ASG-unique pathways found in enrichment results.")
  }

  missing_terms <- dplyr::bind_rows(lapply(names(selected_terms), function(tp) {
    observed <- dot_df$term.id[dot_df$timepoint == tp]
    missing <- setdiff(selected_terms[[tp]], observed)
    if (length(missing) == 0) {
      return(NULL)
    }
    data.frame(timepoint = tp, missing_term_id = missing)
  }))

  if (nrow(missing_terms) > 0) {
    warning(
      "Some selected terms were not plotted because they were absent from ASG-unique enrichment results:\n",
      paste(paste(missing_terms$timepoint, missing_terms$missing_term_id, sep = ": "), collapse = "\n")
    )
  }

  dot_df <- dot_df |>
    dplyr::group_by(timepoint) |>
    dplyr::arrange(fdr, dplyr::desc(number_in_list), .by_group = TRUE) |>
    dplyr::mutate(
      term_order = dplyr::row_number()
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      timepoint = factor(timepoint, levels = names(selected_terms)),
      term_plot = paste(timepoint, term_clean, sep = "___")
    )

  term_levels <- dot_df |>
    dplyr::arrange(timepoint, dplyr::desc(term_order)) |>
    dplyr::pull(term_plot)

  tx_strip_cols <- .get_timepoint_shaded_colors(
    timepoints = names(selected_terms),
    labels = c("Transcript Level", "Overlap", "Gene Level")
  )

  strip_fills <- vapply(
    names(selected_terms),
    function(tp) unname(tx_strip_cols[[tp]][["Transcript Level"]]),
    character(1)
  )

  fig5_d <- ggplot2::ggplot(
    dot_df,
    ggplot2::aes(
      x = gene_ratio,
      y = factor(term_plot, levels = term_levels)
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(size = number_in_list, color = fdr)
    ) +
    ggh4x::facet_wrap2(
      ~timepoint,
      nrow = 1,
      scales = "free_y",
      labeller = ggplot2::labeller(
        timepoint = c(
          H1 = "Hypoxia (1H)",
          H3 = "Hypoxia (3H)",
          H24 = "Hypoxia (24H)"
        )
      ),
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = unname(strip_fills),
          color = unname(grey_palette_colors[["black"]])
        ),
        text_x = ggh4x::elem_list_text(
          color = unname(grey_palette_colors[["black"]]),
          face = "bold"
        )
      )
    ) +
    ggplot2::scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
    ggplot2::scale_x_continuous(labels = scales::percent) +
    ggplot2::scale_color_gradient(
      low = fgsea_colors_up[3],
      high = grey_palette_colors[["mid"]],
      labels = scales::scientific
    ) +
    ggplot2::labs(
      x = "Gene Ratio",
      y = NULL,
      size = "Genes",
      color = "FDR"
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(
        color = unname(grey_palette_colors["light"]),
        linewidth = 0.5
      ),
      panel.grid.major.x = ggplot2::element_line(
        color = unname(grey_palette_colors["light"]),
        linewidth = 0.5
      ),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = base_size - 3.5),
      axis.text.x = ggplot2::element_text(size = base_size - 3),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title.align = 0.5,
      legend.title = ggplot2::element_text(size = base_size),
      legend.text = ggplot2::element_text(size = base_size - 3.5),
      legend.spacing.x = grid::unit(0.5, "cm"),
      legend.box.spacing = grid::unit(0.35, "cm")
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        nrow = 1
      ),
      color = ggplot2::guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        direction = "horizontal",
        barwidth = grid::unit(5, "cm"),
        barheight = grid::unit(0.5, "cm"),
        label.theme = ggplot2::element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = base_size - 3
        )
      )
    )

  ggplot2::ggsave(
    filename = file.path(outdir, filename),
    plot = fig5_d,
    width = 18,
    height = 9
  )

  fig5_d
}

plot_fig_5_all <- function(
  dexseq_results,
  drimseq_results,
  suppa_results,
  deseq_results,
  representation_res = NULL,
  enrichment_res = NULL,
  outdir = NULL,
  use_cached_results = TRUE,
  regenerate_analyses = FALSE,
  regenerate_plots = TRUE,
  save_tables = FALSE,
  enrichment_cutoff = 0.10,
  annot_dataset = "GO:0008150",
  force_enrich = FALSE,
  force_mapping = FALSE,
  force_obo = FALSE,
  force_term2gene = FALSE
) {
  .check_pkg(c("ggplot2", "eulerr", "ggh4x"))

  if (is.null(outdir)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_5")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  representation_rds <- file.path(
    outdir,
    "fig_5_splicing_pathway_representation_analysis.rds"
  )

  enrichment_rds <- file.path(
    outdir,
    "fig_5_splicing_pathway_enrichment_analysis.rds"
  )

  if (is.null(representation_res)) {
    if (
      isTRUE(use_cached_results) &&
        !isTRUE(regenerate_analyses) &&
        file.exists(representation_rds)
    ) {
      message("Loading cached representation analysis: ", representation_rds)
      representation_res <- readRDS(representation_rds)
    } else {
      message("Running pathway representation analysis...")
      representation_res <- run_splicing_pathway_representation_analysis(
        dexseq_results = dexseq_results,
        drimseq_results = drimseq_results,
        suppa_results = suppa_results,
        deseq_results = deseq_results,
        outdir = outdir,
        save_tables = save_tables,
        force_mapping = force_mapping,
        force_obo = force_obo,
        force_term2gene = force_term2gene
      )
    }
  }

  if (is.null(enrichment_res)) {
    if (
      isTRUE(use_cached_results) &&
        !isTRUE(regenerate_analyses) &&
        file.exists(enrichment_rds)
    ) {
      message("Loading cached pathway enrichment analysis: ", enrichment_rds)
      enrichment_res <- readRDS(enrichment_rds)
    } else {
      message("Running transcript-level pathway enrichment analysis...")
      enrichment_res <- run_splicing_pathway_enrichment_analysis(
        dexseq_results = dexseq_results,
        drimseq_results = drimseq_results,
        suppa_results = suppa_results,
        deseq_results = deseq_results,
        enrichment_cutoff = enrichment_cutoff,
        annot_dataset = annot_dataset,
        outdir = outdir,
        save_tables = save_tables,
        force_enrich = force_enrich
      )
    }
  }

  figs <- NULL

  if (isTRUE(regenerate_plots)) {
    message("Generating Figure 5 panels...")

    fig5_a <- plot_clustered_pathway_representation(
      representation_res = representation_res,
      outdir = outdir
    )

    fig5_b <- plot_splicing_pathway_overlap_across_timepoints(
      res = enrichment_res,
      outdir = outdir
    )

    fig5_c <- plot_splicing_vs_deg_pathway_overlap(
      res = enrichment_res,
      outdir = outdir
    )

    fig5_d <- plot_transcript_unique_pathway_dotplot_curated(
      res = enrichment_res,
      outdir = outdir
    )

    figs <- list(
      A = fig5_a,
      B = fig5_b,
      C = fig5_c,
      D = fig5_d
    )
  }

  invisible(list(
    representation_res = representation_res,
    enrichment_res = enrichment_res,
    figs = figs,
    paths = list(
      outdir = outdir,
      representation_rds = representation_rds,
      enrichment_rds = enrichment_rds,
      fig5_a = file.path(outdir, "fig_5_A_clustered_pathway_representation.pdf"),
      fig5_b = file.path(outdir, "fig_5_B_splicing_pathway_overlap_across_timepoints.pdf"),
      fig5_c = file.path(outdir, "fig_5_C_splicing_vs_deg_pathway_overlap.pdf"),
      fig5_d = file.path(outdir, "fig_5_D_transcript_unique_pathways_curated_dotplot.pdf")
    )
  ))
}
