# ============================================================
# plot_sup_fig_5.R
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

.reference_conditions <- function() {
  c("C1", "H1", "H3", "H24")
}

.unique_pair_ids <- function(all_conditions = .reference_conditions()) {
  utils::combn(
    all_conditions,
    2,
    FUN = function(x) paste(x, collapse = "_vs_")
  )
}

.make_pair_id <- function(a, b, all_conditions = .reference_conditions()) {
  ord <- match(c(a, b), all_conditions)
  if (anyNA(ord)) {
    stop("Unknown condition in comparison: ", paste(c(a, b), collapse = ", "))
  }
  paste(c(a, b)[order(ord)], collapse = "_vs_")
}

.get_unique_test_levels <- function(all_conditions = .reference_conditions(), ref_level) {
  ref_i <- match(ref_level, all_conditions)

  if (is.na(ref_i) || ref_i >= length(all_conditions)) {
    return(character(0))
  }

  all_conditions[(ref_i + 1):length(all_conditions)]
}

.comparison_label <- function(pair_id) {
  gsub("_vs_", " -> ", pair_id, fixed = TRUE)
}

.contrast_group_from_pair_id <- function(pair_id) {
  ifelse(grepl("(^C1_vs_|_vs_C1$)", pair_id), "Normoxia vs Hypoxia", "Hypoxia vs Hypoxia")
}

.require_columns <- function(df, cols, df_name = deparse(substitute(df))) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(
      df_name,
      " is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }
  invisible(TRUE)
}

.standardize_gene_column <- function(df) {
  if ("ensgene" %in% names(df)) {
    return(df)
  }
  if ("gene_id" %in% names(df)) {
    df$ensgene <- df$gene_id
    return(df)
  }
  if ("gene" %in% names(df)) {
    df$ensgene <- df$gene
    return(df)
  }
  stop("Could not find a gene identifier column. Expected one of: ensgene, gene_id, gene.")
}

.standardize_reflevel_results <- function(
  ref_results,
  tool,
  all_conditions = .reference_conditions(),
  sig_slot = NULL
) {
  .check_pkg(c("dplyr", "purrr"))

  if (is.null(ref_results) || length(ref_results) == 0) {
    stop("ref_results is empty for ", tool, ".")
  }

  wanted_pair_ids <- .unique_pair_ids(all_conditions)

  out <- purrr::imap_dfr(ref_results, function(x, ref_name) {
    ref_level <- NULL

    if (!is.null(x$meta$ref_level)) {
      ref_level <- x$meta$ref_level
    } else {
      ref_level <- ref_name
    }

    if (!is.null(sig_slot)) {
      for (slot in sig_slot) {
        if (!is.null(x[[slot]])) {
          x <- x[[slot]]
          break
        }
      }
    }

    if (is.list(x) && !is.data.frame(x)) {
      if (!is.null(x$results$gene_sig_all)) {
        df <- x$results$gene_sig_all
      } else if (!is.null(x$gene_sig_all)) {
        df <- x$gene_sig_all
      } else if (!is.null(x$results)) {
        df <- x$results
      } else {
        stop("Could not identify a results data frame for ", tool, ", ref = ", ref_level)
      }
    } else {
      df <- x
    }

    df <- as.data.frame(df, stringsAsFactors = FALSE)
    df <- .standardize_gene_column(df)

    if (!("timepoint" %in% names(df))) {
      if ("test_level" %in% names(df)) {
        df$timepoint <- df$test_level
      } else if ("condition" %in% names(df)) {
        df$timepoint <- df$condition
      } else {
        stop("Could not find comparison/test condition column for ", tool, ", ref = ", ref_level)
      }
    }

    df$ref_level <- ref_level
    df$tool <- tool
    df
  })

  out |>
    dplyr::mutate(
      pair_id = purrr::map2_chr(
        ref_level,
        timepoint,
        ~ .make_pair_id(.x, .y, all_conditions = all_conditions)
      ),
      comparison = .comparison_label(pair_id),
      contrast_group = .contrast_group_from_pair_id(pair_id)
    ) |>
    dplyr::filter(pair_id %in% wanted_pair_ids) |>
    dplyr::distinct(tool, pair_id, comparison, contrast_group, ensgene, .keep_all = TRUE)
}

.make_upset_input <- function(pairwise_df) {
  .check_pkg("dplyr")

  .require_columns(pairwise_df, c("pair_id", "ensgene"), "pairwise_df")

  split(pairwise_df$ensgene, pairwise_df$pair_id) |>
    lapply(unique)
}

# --------------------------------------------------
# Reference-level Analyses
# --------------------------------------------------

run_deseq_all_unique_reflevels <- function(
  all_conditions = .reference_conditions(),
  force = FALSE,
  ...
) {
  out <- list()

  for (ref_level in all_conditions) {
    test_levels <- .get_unique_test_levels(all_conditions, ref_level)

    if (length(test_levels) == 0) {
      message("[DESeq2] No new unique comparisons for reference level = ", ref_level)
      next
    }

    outdir <- file.path(.get_results_dir(), "analysis", "deseq")
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    out_rds <- file.path(outdir, paste0("deseq_results.", ref_level, ".rds"))
    out_xlsx <- file.path(outdir, paste0("deseq_results.", ref_level, ".xlsx"))

    if (file.exists(out_rds) && !force) {
      message("[DESeq2] Cached RDS exists; loading reference level = ", ref_level)
      out[[ref_level]] <- readRDS(out_rds)
      next
    }

    message(
      "[DESeq2] Running reference level = ", ref_level,
      " comparisons: ",
      paste(test_levels, "vs", ref_level, collapse = ", ")
    )

    out[[ref_level]] <- run_deseq(
      ref_level = ref_level,
      timepoints = test_levels,
      force = force,
      out_rds = out_rds,
      out_xlsx = out_xlsx,
      ...
    )
  }

  invisible(out)
}

run_dexseq_all_unique_reflevels <- function(
  all_conditions = .reference_conditions(),
  mode_tag = "summarizeOverlaps.multiOverlap.all",
  force = FALSE,
  ...
) {
  out <- list()

  for (ref_level in all_conditions) {
    test_levels <- .get_unique_test_levels(all_conditions, ref_level)

    if (length(test_levels) == 0) {
      message("[DEXSeq] No new unique comparisons for reference level = ", ref_level)
      next
    }

    outdir <- file.path(.get_results_dir(), "analysis", "dexseq")
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    out_rds <- file.path(outdir, paste0("dexseq_results.", ref_level, ".", mode_tag, ".rds"))
    out_xlsx <- file.path(outdir, paste0("dexseq_results.", ref_level, ".", mode_tag, ".xlsx"))

    if (file.exists(out_rds) && !force) {
      message("[DEXSeq] Loading cached default-mode RDS for reference level = ", ref_level, ": ", basename(out_rds))
      out[[ref_level]] <- readRDS(out_rds)
      next
    }

    message(
      "[DEXSeq] Running reference level = ", ref_level,
      " comparisons: ",
      paste(test_levels, "vs", ref_level, collapse = ", ")
    )

    out[[ref_level]] <- run_dexseq(
      ref_level = ref_level,
      timepoints = test_levels,
      force = force,
      out_rds = out_rds,
      out_xlsx = out_xlsx,
      ...
    )

    gc()
  }

  invisible(out)
}

run_drimseq_all_unique_reflevels <- function(
  all_conditions = .reference_conditions(),
  force = FALSE,
  ...
) {
  out <- list()

  for (ref_level in all_conditions) {
    test_levels <- .get_unique_test_levels(all_conditions, ref_level)

    if (length(test_levels) == 0) {
      message("[DRIMSeq] No new unique comparisons for reference level = ", ref_level)
      next
    }

    outdir <- file.path(.get_results_dir(), "analysis", "drimseq")
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    out_rds <- file.path(outdir, paste0("drimseq_results.", ref_level, ".rds"))
    out_xlsx <- file.path(outdir, paste0("drimseq_results.", ref_level, ".xlsx"))

    if (file.exists(out_rds) && !force) {
      message("[DRIMSeq] Cached RDS exists; loading reference level = ", ref_level)
      out[[ref_level]] <- readRDS(out_rds)
      next
    }

    message(
      "[DRIMSeq] Running reference level = ", ref_level,
      " comparisons: ",
      paste(test_levels, "vs", ref_level, collapse = ", ")
    )

    out[[ref_level]] <- run_drimseq(
      ref_level = ref_level,
      timepoints = test_levels,
      force = force,
      out_rds = out_rds,
      out_xlsx = out_xlsx,
      ...
    )

    gc()
  }

  invisible(out)
}

# --------------------------------------------------
# Main Analysis Object Builder
# --------------------------------------------------

make_reference_level_plot_data <- function(
  deseq_ref_results = NULL,
  dexseq_ref_results = NULL,
  drimseq_ref_results = NULL,
  suppa_ref_results = NULL,
  all_conditions = .reference_conditions()
) {
  .check_pkg(c("dplyr", "purrr"))

  out <- list()

  if (!is.null(deseq_ref_results)) {
    out$deseq2 <- .standardize_reflevel_results(
      ref_results = deseq_ref_results,
      tool = "DESeq2",
      all_conditions = all_conditions
    )
  }

  if (!is.null(dexseq_ref_results)) {
    out$dexseq <- .standardize_reflevel_results(
      ref_results = dexseq_ref_results,
      tool = "DEXSeq",
      all_conditions = all_conditions
    )
  }

  if (!is.null(drimseq_ref_results)) {
    out$drimseq <- .standardize_reflevel_results(
      ref_results = drimseq_ref_results,
      tool = "DRIMSeq",
      all_conditions = all_conditions
    )
  }

  if (!is.null(suppa_ref_results)) {
    out$suppa2 <- .standardize_reflevel_results(
      ref_results = suppa_ref_results,
      tool = "SUPPA2",
      all_conditions = all_conditions
    )
  }

  if (length(out) == 0) {
    stop("No reference-level result objects were provided.")
  }

  pairwise <- dplyr::bind_rows(out)

  invisible(list(
    pairwise = pairwise,
    by_tool = out,
    all_conditions = all_conditions
  ))
}

compute_reference_dependency <- function(pairwise_df) {
  .check_pkg("dplyr")

  pairwise_df |>
    dplyr::distinct(ensgene, pair_id) |>
    dplyr::mutate(is_normoxia = grepl("(^C1_vs_|_vs_C1$)", pair_id)) |>
    dplyr::group_by(ensgene) |>
    dplyr::summarize(
      any_normoxia = any(is_normoxia),
      any_hypoxia = any(!is_normoxia),
      category = dplyr::case_when(
        any_normoxia & any_hypoxia ~ "Shared",
        any_normoxia & !any_hypoxia ~ "Normoxia\nContrast Only",
        !any_normoxia & any_hypoxia ~ "Hypoxia\nContrast Only",
        TRUE ~ NA_character_
      ),
      .groups = "drop"
    ) |>
    dplyr::count(category, name = "n_genes") |>
    dplyr::mutate(prop = n_genes / sum(n_genes))
}

# --------------------------------------------------
# Plotting Helpers
# --------------------------------------------------

.tool_color <- function(tool_name, shade = "base") {
  base_cols <- splicing_tool_colors

  if (!("DESeq2" %in% names(base_cols))) {
    base_cols <- c("DESeq2" = unname(category_base_colors[["Category 2"]]), base_cols)
  }

  if (!(tool_name %in% names(base_cols))) {
    base_cols <- c(base_cols, "Fallback" = unname(grey_palette_colors[["dark"]]))
    tool_name <- "Fallback"
  }

  shaded <- .generate_shaded_palette(
    base_colors = base_cols,
    labels = c("light", "base", "dark")
  )

  unname(shaded[[tool_name]][[shade]])
}

.reference_level_labels <- function() {
  c(
    C1 = "Normoxia",
    H1 = "Hypoxia (1H)",
    H3 = "Hypoxia (3H)",
    H24 = "Hypoxia (24H)"
  )
}

.reference_level_colors <- function(shade = "base") {
  ref_labels <- .reference_level_labels()
  base_cols <- timepoint_base_colors[unname(ref_labels)]
  names(base_cols) <- names(ref_labels)

  shaded <- .generate_shaded_palette(
    base_colors = base_cols,
    labels = c("light", "base", "dark")
  )

  stats::setNames(
    vapply(shaded, function(x) unname(x[[shade]]), character(1)),
    names(shaded)
  )
}

.dependency_order <- function() {
  c(
    "Normoxia\nContrast Only",
    "Shared",
    "Hypoxia\nContrast Only"
  )
}

.tool_dependency_colors <- function(tool_name) {
  stats::setNames(
    c(
      .tool_color(tool_name, shade = "dark"),
      .tool_color(tool_name, shade = "base"),
      .tool_color(tool_name, shade = "light")
    ),
    .dependency_order()
  )
}

.pair_dependency_category <- function(pair_id) {
  ifelse(
    grepl("(^C1_vs_|_vs_C1$)", pair_id),
    "Normoxia\nContrast Only",
    "Hypoxia\nContrast Only"
  )
}

.dependency_legend_labels <- function() {
  c(
    "Normoxia\nContrast Only" = "Normoxia Contrast Only",
    "Shared" = "Shared",
    "Hypoxia\nContrast Only" = "Hypoxia Contrast Only"
  )
}

.make_set_metadata <- function(set_order) {
  data.frame(
    group = set_order,
    ref_level = sub("_vs_.*$", "", set_order),
    stringsAsFactors = FALSE
  )
}

# --------------------------------------------------
# Plotting Functions
# --------------------------------------------------

plot_reference_level_upset <- function(
  pairwise_df,
  tool_name = NULL,
  outdir = NULL,
  filename = NULL,
  title = NULL,
  width = 7,
  height = 5,
  min_size = 10,
  base_size = 12,
  save = TRUE
) {
  .check_pkg(c("ComplexUpset", "ggplot2", "dplyr"))

  outdir <- .set_outdir(outdir, subdir = "sup_fig_5")

  if (is.null(tool_name) && "tool" %in% names(pairwise_df)) {
    tool_name <- unique(pairwise_df$tool)
    if (length(tool_name) != 1) tool_name <- NULL
  }

  if (is.null(filename)) {
    filename <- paste0(
      "reference_level_upset_",
      tolower(gsub("[^A-Za-z0-9]+", "_", tool_name)),
      ".pdf"
    )
  }

  tool_col <- .tool_color(tool_name)

  inactive_dot_col <- .make_shades_one(
    base_color = tool_col,
    labels = c("light", "base", "dark")
  )[["light"]]

  dep_cols <- .tool_dependency_colors(tool_name)

  sets <- .make_upset_input(pairwise_df)
  all_genes <- sort(unique(unlist(sets, use.names = FALSE)))

  set_order <- .unique_pair_ids(.reference_conditions())
  set_order <- intersect(set_order, names(sets))

  set_labels <- stats::setNames(
    .comparison_label(set_order),
    set_order
  )

  display_order <- unname(set_labels[set_order])

  upset_df <- data.frame(ensgene = all_genes, stringsAsFactors = FALSE)
  for (set_name in set_order) {
    display_name <- set_labels[[set_name]]
    upset_df[[display_name]] <- upset_df$ensgene %in% sets[[set_name]]
  }

  set_metadata <- data.frame(
    set = rev(display_order),
    dependency = .pair_dependency_category(rev(set_order)),
    stringsAsFactors = FALSE
  )

  set_metadata$dependency <- factor(
    set_metadata$dependency,
    levels = .dependency_order()
  )

  fig <- suppressWarnings(ComplexUpset::upset(
    upset_df,
    intersect = rev(display_order),
    name = "Contrast Overlap",
    min_size = min_size,
    sort_sets = FALSE,
    sort_intersections_by = "cardinality",
    height_ratio = 0.9,
    set_sizes = FALSE,
    base_annotations = list(
      "Significant Genes" =
        ComplexUpset::intersection_size(
          mapping = ggplot2::aes(fill = "bars_color"),
          text = list(vjust = -0.35, hjust = -0.05, size = 3, angle = 30),
          counts = TRUE
        ) +
          ggplot2::scale_fill_manual(
            guide = "none",
            values = c("bars_color" = tool_col)
          )
    ),
    matrix = ComplexUpset::intersection_matrix(
      geom = ggplot2::geom_point(size = 3),
      segment = ggplot2::geom_segment(linewidth = 1.25),
      outline_color = list(
        active = tool_col,
        inactive = inactive_dot_col
      )
    ),
    stripes = ComplexUpset::upset_stripes(
      mapping = ggplot2::aes(color = dependency),
      data = set_metadata,
      colors = dep_cols
    ),
    themes = ComplexUpset::upset_modify_themes(
      list(
        "Significant Genes" = ggplot2::theme(
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_text(size = base_size),
          panel.grid.major = ggplot2::element_line(color = unname(grey_palette_colors["light"])),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_rect(fill = "white", colour = NA),
          plot.background = ggplot2::element_rect(fill = "white", colour = NA)
        )
      )
    )
  ) +
    ggplot2::labs(title = title) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 18, face = "bold"),
      text = ggplot2::element_text(size = base_size)
    ))

  if (isTRUE(save)) {
    ggplot2::ggsave(
      filename = file.path(outdir, filename),
      plot = fig,
      width = width,
      height = height
    )
  }

  invisible(fig)
}

plot_reference_set_sizes <- function(
  pairwise_df,
  tool_name = NULL,
  outdir = NULL,
  filename = NULL,
  title = NULL,
  base_size = 12,
  width = 5,
  height = 5,
  save = TRUE
) {
  .check_pkg(c("dplyr", "ggplot2", "scales", "ggh4x"))

  outdir <- .set_outdir(outdir, subdir = "sup_fig_5")

  if (is.null(tool_name) && "tool" %in% names(pairwise_df)) {
    tool_name <- unique(pairwise_df$tool)
    if (length(tool_name) != 1) tool_name <- NULL
  }

  if (is.null(filename)) {
    filename <- paste0(
      "reference_level_set_sizes_",
      tolower(gsub("[^A-Za-z0-9]+", "_", tool_name)),
      ".pdf"
    )
  }

  category_order <- c(
    "Normoxia\nContrast Only",
    "Hypoxia\nContrast Only"
  )

  pair_order <- .unique_pair_ids(.reference_conditions())
  pair_labels <- stats::setNames(
    paste0(.comparison_label(pair_order), "\n\n"),
    pair_order
  )

  tool_cols <- .tool_dependency_colors(tool_name)

  df <- pairwise_df |>
    dplyr::distinct(pair_id, comparison, contrast_group, ensgene) |>
    dplyr::count(pair_id, comparison, contrast_group, name = "n_genes") |>
    dplyr::mutate(
      category = .pair_dependency_category(pair_id),
      category = factor(category, levels = category_order),
      pair_id = factor(pair_id, levels = pair_order)
    )

  strip_cols <- tool_cols[category_order]

  fig <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = pair_id, y = n_genes, fill = category)
  ) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::comma(n_genes)),
      vjust = -0.35,
      size = 3
    ) +
    ggh4x::facet_grid2(
      ~category,
      scales = "free_x",
      space = "free_x",
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = unname(strip_cols),
          color = unname(grey_palette_colors[["black"]])
        ),
        text_x = ggh4x::elem_list_text(
          color = "white",
          face = "bold"
        )
      )
    ) +
    ggplot2::scale_x_discrete(
      drop = TRUE,
      labels = pair_labels
    ) +
    ggplot2::scale_fill_manual(values = tool_cols, drop = FALSE) +
    ggplot2::scale_y_continuous(
      labels = scales::comma,
      expand = ggplot2::expansion(mult = c(0, 0.15))
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Significant Genes",
      title = title
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold")
    )

  if (isTRUE(save)) {
    ggplot2::ggsave(
      filename = file.path(outdir, filename),
      plot = fig,
      width = width,
      height = height
    )
  }

  invisible(fig)
}

plot_reference_dependency <- function(
  pairwise_df,
  tool_name = NULL,
  outdir = NULL,
  filename = NULL,
  title = NULL,
  base_size = 12,
  width = 3.5,
  height = 5,
  save = TRUE
) {
  .check_pkg(c("dplyr", "ggplot2", "scales", "tidyr", "ggh4x"))

  outdir <- .set_outdir(outdir, subdir = "sup_fig_5")

  if (is.null(tool_name) && "tool" %in% names(pairwise_df)) {
    tool_name <- unique(pairwise_df$tool)
    if (length(tool_name) != 1) tool_name <- NULL
  }

  if (is.null(filename)) {
    filename <- paste0(
      "reference_dependency_",
      tolower(gsub("[^A-Za-z0-9]+", "_", tool_name)),
      ".pdf"
    )
  }

  category_order <- .dependency_order()
  dep_cols <- .tool_dependency_colors(tool_name)

  dep_df <- compute_reference_dependency(pairwise_df) |>
    dplyr::mutate(category = factor(category, levels = category_order)) |>
    tidyr::complete(category = factor(category_order, levels = category_order), fill = list(n_genes = 0, prop = 0))

  dep_df$panel_label <- "Contrast\nDependency"
  tool_col <- .tool_color(tool_name)

  fig <- ggplot2::ggplot(dep_df, ggplot2::aes(x = category, y = prop, fill = category)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::percent(prop, accuracy = 1)),
      vjust = -0.35,
      size = 3
    ) +
    ggh4x::facet_grid2(
      ~panel_label,
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = tool_col,
          color = unname(grey_palette_colors["black"])
        ),
        text_x = ggh4x::elem_list_text(
          color = unname(grey_palette_colors["light"]),
          face = "bold"
        )
      )
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(),
      expand = ggplot2::expansion(mult = c(0, 0.15))
    ) +
    ggplot2::scale_fill_manual(values = dep_cols, drop = FALSE) +
    ggplot2::labs(
      x = NULL,
      y = "Significant Genes (%)",
      title = title
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold")
    )

  if (isTRUE(save)) {
    ggplot2::ggsave(
      filename = file.path(outdir, filename),
      plot = fig,
      width = width,
      height = height
    )
  }

  invisible(fig)
}

plot_dependency_top_legend_one <- function(
  tool_name,
  base_size = 15,
  title_size = 18,
  key_size_cm = 0.75
) {
  .check_pkg(c("ggplot2", "cowplot"))

  cats <- .dependency_order()
  cols <- .tool_dependency_colors(tool_name)

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
      labels = .dependency_legend_labels(),
      name = tool_name,
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
      legend.title = ggplot2::element_text(size = title_size, hjust = 0.5),
      legend.text = ggplot2::element_text(size = base_size),
      legend.key.width = grid::unit(key_size_cm, "cm"),
      legend.key.height = grid::unit(key_size_cm, "cm"),
      legend.spacing.x = grid::unit(0.25, "cm"),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )

  cowplot::get_legend(legend_plot)
}

plot_dependency_top_legend <- function(
  base_size = 15,
  title_size = 18,
  key_size_cm = 0.75
) {
  .check_pkg("cowplot")

  deseq_legend <- plot_dependency_top_legend_one(
    tool_name = "DESeq2",
    base_size = base_size,
    title_size = title_size,
    key_size_cm = key_size_cm
  )

  dexseq_legend <- plot_dependency_top_legend_one(
    tool_name = "DEXSeq",
    base_size = base_size,
    title_size = title_size,
    key_size_cm = key_size_cm
  )

  cowplot::plot_grid(
    deseq_legend,
    dexseq_legend,
    nrow = 1,
    align = "h",
    rel_widths = c(1, 1)
  )
}

plot_dependency_bottom_legend_one <- function(
  tool_name,
  base_size = 15,
  title_size = 18,
  key_size_cm = 0.75
) {
  .check_pkg(c("ggplot2", "cowplot"))

  cats <- .dependency_order()
  cols <- .tool_dependency_colors(tool_name)

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
      labels = .dependency_legend_labels(),
      name = tool_name,
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
      legend.title = ggplot2::element_text(size = title_size, hjust = 0.5),
      legend.text = ggplot2::element_text(size = base_size),
      legend.key.width = grid::unit(key_size_cm, "cm"),
      legend.key.height = grid::unit(key_size_cm, "cm"),
      legend.spacing.x = grid::unit(0.25, "cm"),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )

  cowplot::get_legend(legend_plot)
}

plot_dependency_bottom_legend <- function(
  base_size = 15,
  title_size = 18,
  key_size_cm = 0.75
) {
  .check_pkg("cowplot")

  drimseq_legend <- plot_dependency_bottom_legend_one(
    tool_name = "DRIMSeq",
    base_size = base_size,
    title_size = title_size,
    key_size_cm = key_size_cm
  )

  suppa2_legend <- plot_dependency_bottom_legend_one(
    tool_name = "SUPPA2",
    base_size = base_size,
    title_size = title_size,
    key_size_cm = key_size_cm
  )

  cowplot::plot_grid(
    drimseq_legend,
    suppa2_legend,
    nrow = 1,
    align = "h",
    rel_widths = c(1, 1)
  )
}

plot_reference_level_tool_panel <- function(
  pairwise_df,
  tool_name = NULL,
  outdir = NULL,
  filename = NULL,
  title = NULL,
  upset_width = 4.25,
  set_size_width = 2.75,
  dependency_width = 1.5,
  height = 5.5,
  min_size = 15,
  save = FALSE
) {
  .check_pkg(c("patchwork", "ggplot2"))

  outdir <- .set_outdir(outdir, subdir = "sup_fig_5")

  if (is.null(tool_name) && "tool" %in% names(pairwise_df)) {
    tool_name <- unique(pairwise_df$tool)
    if (length(tool_name) != 1) tool_name <- NULL
  }

  if (is.null(title)) {
    title <- tool_name
  }

  if (is.null(filename)) {
    filename <- paste0(
      "reference_level_panel_",
      tolower(gsub("[^A-Za-z0-9]+", "_", tool_name)),
      ".pdf"
    )
  }

  upset_plot <- plot_reference_level_upset(
    pairwise_df = pairwise_df,
    tool_name = tool_name,
    outdir = outdir,
    title = NULL,
    width = upset_width,
    height = height,
    min_size = min_size,
    save = FALSE
  )

  set_size_plot <- plot_reference_set_sizes(
    pairwise_df = pairwise_df,
    tool_name = tool_name,
    outdir = outdir,
    title = NULL,
    width = set_size_width,
    height = height,
    save = FALSE
  )

  dependency_plot <- plot_reference_dependency(
    pairwise_df = pairwise_df,
    tool_name = tool_name,
    outdir = outdir,
    title = NULL,
    width = dependency_width,
    height = height,
    save = FALSE
  )

  panel <- patchwork::wrap_elements(full = upset_plot) +
    patchwork::wrap_elements(full = set_size_plot) +
    patchwork::wrap_elements(full = dependency_plot) +
    patchwork::plot_layout(
      widths = c(upset_width, set_size_width, dependency_width)
    )

  if (isTRUE(save)) {
    ggplot2::ggsave(
      filename = file.path(outdir, filename),
      plot = panel,
      width = upset_width + set_size_width + dependency_width,
      height = height
    )
  }

  invisible(list(
    panel = panel,
    upset = upset_plot,
    set_sizes = set_size_plot,
    dependency = dependency_plot
  ))
}

plot_reference_level_all_tools_panel <- function(
  plot_data,
  outdir = NULL,
  filename = "sup_fig_5_reference_level_analysis.pdf",
  min_size = 10
) {
  .check_pkg(c("patchwork", "ggplot2", "cowplot"))

  outdir <- .set_outdir(outdir, subdir = "sup_fig_5")

  if (is.null(plot_data$by_tool)) {
    stop("Expected plot_data from make_reference_level_plot_data().")
  }

  panels <- lapply(names(plot_data$by_tool), function(tool_key) {
    pairwise_df <- plot_data$by_tool[[tool_key]]
    tool_name <- unique(pairwise_df$tool)
    if (length(tool_name) != 1) tool_name <- tool_key

    res <- plot_reference_level_tool_panel(
      pairwise_df = pairwise_df,
      tool_name = tool_name,
      outdir = outdir,
      title = NULL,
      min_size = min_size
    )

    title_strip <- plot_tool_row_title(tool_name)

    title_strip / patchwork::wrap_elements(full = res$panel) +
      patchwork::plot_layout(heights = c(0.25, 5))
  })

  names(panels) <- names(plot_data$by_tool)

  combined_body <- patchwork::wrap_plots(
    panels,
    ncol = 1
  )

  top_legend <- plot_dependency_top_legend(
    base_size = 15,
    title_size = 18,
    key_size_cm = 0.75
  )

  bottom_legend <- plot_dependency_bottom_legend(
    base_size = 15,
    title_size = 18,
    key_size_cm = 0.75
  )

  legend_panel <- patchwork::wrap_plots(
    patchwork::wrap_elements(full = top_legend),
    patchwork::wrap_elements(full = bottom_legend),
    ncol = 1,
    heights = c(1, 1)
  )

  body_height <- 5.5 * length(panels)
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
    width = 20,
    height = body_height + legend_height
  )

  invisible(list(
    combined = combined,
    panels = panels
  ))
}

# --------------------------------------------------
# Convenience Wrapper (Plotting Only)
# --------------------------------------------------

plot_reference_level_supplement_all <- function(
  plot_data,
  outdir = NULL,
  make_tool_panels = FALSE,
  make_combined_panel = TRUE,
  make_upset = FALSE,
  make_dependency = FALSE,
  min_size = 10
) {
  outdir <- .set_outdir(outdir, subdir = "sup_fig_5")

  if (is.null(plot_data$by_tool)) {
    stop("Expected plot_data from make_reference_level_plot_data().")
  }

  figs <- list()

  if (isTRUE(make_tool_panels)) {
    figs$tool_panels <- list()
  }

  for (tool_key in names(plot_data$by_tool)) {
    pairwise_df <- plot_data$by_tool[[tool_key]]
    tool_name <- unique(pairwise_df$tool)
    if (length(tool_name) != 1) tool_name <- tool_key

    figs[[tool_key]] <- list()

    if (isTRUE(make_tool_panels)) {
      figs$tool_panels[[tool_key]] <- plot_reference_level_tool_panel(
        pairwise_df = pairwise_df,
        tool_name = tool_name,
        outdir = outdir,
        min_size = min_size
      )
    }

    if (isTRUE(make_upset)) {
      figs[[tool_key]]$upset <- plot_reference_level_upset(
        pairwise_df = pairwise_df,
        tool_name = tool_name,
        outdir = outdir,
        min_size = min_size
      )
    }

    if (isTRUE(make_dependency)) {
      figs[[tool_key]]$dependency <- plot_reference_dependency(
        pairwise_df = pairwise_df,
        tool_name = tool_name,
        outdir = outdir
      )
    }
  }

  if (isTRUE(make_combined_panel)) {
    figs$combined_panel <- plot_reference_level_all_tools_panel(
      plot_data = plot_data,
      outdir = outdir,
      min_size = min_size
    )
  }

  invisible(figs)
}

# --------------------------------------------------
# Convenience Wrapper (Full Pipeline)
# --------------------------------------------------

plot_sup_fig_5_all <- function(
  all_conditions = .reference_conditions(),
  force = FALSE,
  outdir = NULL,
  plot_data_rds = NULL,
  use_cached_plot_data = TRUE,
  save_plot_data = TRUE,
  regenerate_plots = TRUE,
  min_size = 15
) {
  outdir <- .set_outdir(outdir, subdir = "sup_fig_5")

  if (is.null(plot_data_rds)) {
    plot_data_rds <- file.path(outdir, "sup_fig_5_reference_level_analysis.rds")
  }

  message("[Sup Fig 5] Starting full reference-level analysis + plotting")

  if (isTRUE(use_cached_plot_data) && file.exists(plot_data_rds) && !force) {
    message("[Sup Fig 5] Loading cached plot data: ", plot_data_rds)
    plot_data <- readRDS(plot_data_rds)

    figs <- NULL
    if (isTRUE(regenerate_plots)) {
      message("\n[Sup Fig 5] === Generating plots from cached plot data ===")
      figs <- plot_reference_level_supplement_all(
        plot_data,
        outdir = outdir,
        make_tool_panels = FALSE,
        make_combined_panel = TRUE,
        make_upset = FALSE,
        make_dependency = FALSE,
        min_size = min_size
      )
    }

    message("\n[Sup Fig 5] Completed successfully using cached plot data")

    return(invisible(list(
      deseq = NULL,
      dexseq = NULL,
      drimseq = NULL,
      suppa = NULL,
      plot_data = plot_data,
      figs = figs,
      paths = list(outdir = outdir, plot_data_rds = plot_data_rds)
    )))
  }

  message("\n[Sup Fig 5] === Running DESeq2 ===")
  deseq_ref_results <- run_deseq_all_unique_reflevels(all_conditions = all_conditions, force = force)

  message("\n[Sup Fig 5] === Running DEXSeq ===")
  dexseq_ref_results <- run_dexseq_all_unique_reflevels(all_conditions = all_conditions, force = force)

  message("\n[Sup Fig 5] === Running DRIMSeq ===")
  drimseq_ref_results <- run_drimseq_all_unique_reflevels(all_conditions = all_conditions, force = force)

  message("\n[Sup Fig 5] === Loading SUPPA2 ===")
  suppa_ref_results <- list()

  for (ref_level in all_conditions) {
    test_levels <- .get_unique_test_levels(all_conditions, ref_level)

    if (length(test_levels) == 0) {
      message("[SUPPA2] No new unique comparisons for reference level = ", ref_level)
      next
    }

    message(
      "[SUPPA2] Loading reference level = ", ref_level,
      " comparisons: ",
      paste(test_levels, "vs", ref_level, collapse = ", ")
    )

    suppa_ref_results[[ref_level]] <- load_suppa(
      ref_level = ref_level,
      timepoints = test_levels,
      force = force
    )
  }

  message("\n[Sup Fig 5] === Building Plot Data ===")
  plot_data <- make_reference_level_plot_data(
    deseq_ref_results = deseq_ref_results,
    dexseq_ref_results = dexseq_ref_results,
    drimseq_ref_results = drimseq_ref_results,
    suppa_ref_results = suppa_ref_results,
    all_conditions = all_conditions
  )

  if (isTRUE(save_plot_data)) {
    saveRDS(plot_data, plot_data_rds)
    message("[Sup Fig 5] Saved plot data: ", plot_data_rds)
  }

  figs <- NULL
  if (isTRUE(regenerate_plots)) {
    message("\n[Sup Fig 5] === Generating Plots ===")
    figs <- plot_reference_level_supplement_all(
      plot_data,
      outdir = outdir,
      make_tool_panels = FALSE,
      make_combined_panel = TRUE,
      make_upset = FALSE,
      make_dependency = FALSE,
      min_size = min_size
    )
  }

  message("\n[Sup Fig 5] Completed successfully...")

  invisible(list(
    deseq = deseq_ref_results,
    dexseq = dexseq_ref_results,
    drimseq = drimseq_ref_results,
    suppa = suppa_ref_results,
    plot_data = plot_data,
    figs = figs,
    paths = list(outdir = outdir, plot_data_rds = plot_data_rds)
  ))
}
