# ============================================================
# plot_sup_fig_1.R
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

.suppa_timepoint_order <- function() {
  c("H1", "H3", "H24")
}

.suppa_timepoint_labels <- function() {
  stats::setNames(
    names(timepoint_base_colors)[grepl("Hypoxia", names(timepoint_base_colors))],
    .suppa_timepoint_order()
  )
}

.suppa_timepoint_colors <- function(shade = "light") {
  timepoints <- .suppa_timepoint_order()

  tp_shades <- .get_timepoint_shaded_colors(
    timepoints = timepoints,
    labels = c("light", "base", "dark")
  )

  stats::setNames(
    vapply(
      tp_shades,
      function(x) unname(x[[shade]]),
      character(1)
    ),
    timepoints
  )
}

.suppa_timepoint_strip <- function(shade = "light") {
  .check_pkg("ggh4x")

  ggh4x::strip_themed(
    background_x = ggh4x::elem_list_rect(
      fill = unname(.suppa_timepoint_colors(shade = shade)),
      color = unname(grey_palette_colors["black"])
    ),
    text_x = ggh4x::elem_list_text(
      color = unname(grey_palette_colors["black"]),
      face = "bold"
    )
  )
}

.bind_suppa_sig_events <- function(suppa_results) {
  .check_pkg("dplyr")

  if (!("by_tp" %in% names(suppa_results))) {
    stop("Expected suppa_results$by_tp.")
  }

  timepoint_order <- .suppa_timepoint_order()

  dplyr::bind_rows(
    lapply(timepoint_order, function(tp) {
      x <- suppa_results$by_tp[[tp]]$sig
      if (is.null(x)) {
        return(NULL)
      }

      x <- as.data.frame(x, stringsAsFactors = FALSE)

      if (!("timepoint" %in% names(x))) {
        x$timepoint <- tp
      } else {
        x$timepoint[is.na(x$timepoint) | x$timepoint == ""] <- tp
      }

      x
    })
  ) |>
    dplyr::mutate(
      timepoint = factor(timepoint, levels = timepoint_order),
      symbol = dplyr::if_else(is.na(symbol) | symbol == "", ensgene, symbol),
      abs_dpsi = abs(dpsi)
    )
}

.bind_suppa_events_for_plots <- function(suppa_results, sig_only = FALSE, p_cutoff = 0.05) {
  .check_pkg("dplyr")

  timepoint_order <- .suppa_timepoint_order()

  dplyr::bind_rows(lapply(timepoint_order, function(tp) {
    x <- if (isTRUE(sig_only)) {
      suppa_results$by_tp[[tp]]$sig
    } else if (!is.null(suppa_results$by_tp[[tp]]$full)) {
      suppa_results$by_tp[[tp]]$full
    } else {
      suppa_results$by_tp[[tp]]$sig
    }

    if (is.null(x)) {
      return(NULL)
    }
    x <- as.data.frame(x, stringsAsFactors = FALSE)

    if (!("timepoint" %in% names(x))) x$timepoint <- tp
    x$timepoint[is.na(x$timepoint) | x$timepoint == ""] <- tp

    x
  })) |>
    dplyr::mutate(
      timepoint = factor(timepoint, levels = timepoint_order),
      symbol = dplyr::if_else(is.na(symbol) | symbol == "", ensgene, symbol),
      neglog10_pvalue = -log10(pvalue),
      abs_dpsi = abs(dpsi),
      sig_class = dplyr::case_when(
        pvalue < p_cutoff & abs_dpsi >= 0.30 ~ "Significant, |Delta PSI| >= 0.30",
        pvalue < p_cutoff ~ "Significant",
        abs_dpsi >= 0.30 ~ "|Delta PSI| >= 0.30",
        TRUE ~ "Not significant"
      )
    )
}

.make_suppa_multiplicity_magnitude_df <- function(suppa_results, p_cutoff = 0.05) {
  .check_pkg("dplyr")

  .bind_suppa_sig_events(suppa_results) |>
    dplyr::filter(
      !is.na(ensgene),
      !is.na(event_id),
      is.finite(abs_dpsi),
      !is.na(pvalue),
      pvalue <= p_cutoff
    ) |>
    dplyr::group_by(timepoint, ensgene, symbol) |>
    dplyr::summarize(
      n_sig_events = dplyr::n_distinct(event_id),
      n_event_types = dplyr::n_distinct(event_type),
      max_abs_dpsi = max(abs_dpsi, na.rm = TRUE),
      mean_abs_dpsi = mean(abs_dpsi, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(tool = "SUPPA2")
}

# --------------------------------------------------
# Supplementary Figure 1A - SUPPA2 Gene Overlap
# --------------------------------------------------

plot_suppa_overlap <- function(
  suppa_results,
  outdir = NULL,
  title = "SUPPA2: Gene Overlap"
) {
  .check_pkg("grid")

  if (is.null(outdir)) {
    outdir <- file.path(.get_results_dir(), "plots", "sup_fig_1")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  features_by_tp <- lapply(c("H1", "H3", "H24"), function(tp) {
    if (!is.null(suppa_results$results$gene_sig_by_tp[[tp]])) {
      suppa_results$results$gene_sig_by_tp[[tp]]$ensgene
    } else {
      suppa_results$by_tp[[tp]]$sig$ensgene
    }
  })
  names(features_by_tp) <- c("H1", "H3", "H24")

  fig <- plot_feature_overlap_across_timepoints(
    features_by_tp = features_by_tp,
    title = title,
    shade_label = "TX"
  )

  grid::grid.newpage()
  grid::grid.draw(fig$plot)

  grDevices::pdf(
    file.path(outdir, "sup_fig_1_A_suppa_gene_overlap.pdf"),
    width = 6.5,
    height = 6.5
  )
  grid::grid.newpage()
  grid::grid.draw(fig$plot)
  invisible(grDevices::dev.off())

  invisible(fig)
}

# --------------------------------------------------
# Supplementary Figure 1B - SUPPA2 Multiplicity vs Magnitude
# --------------------------------------------------

plot_suppa_multiplicity_vs_magnitude <- function(
  suppa_results,
  outdir = NULL,
  title = "SUPPA2: Multiplicity and Magnitude",
  base_size = 15,
  p_cutoff = 0.10,
  max_events_to_show = 10,
  label_top_n_by_events = 5,
  label_min_events = 5,
  point_alpha = 0.75,
  point_size = 3
) {
  .check_pkg(c(
    "dplyr",
    "ggplot2",
    "ggrepel",
    "ggh4x",
    "scales"
  ))

  outdir <- .set_outdir(outdir, subdir = "sup_fig_1")
  timepoint_labels <- .suppa_timepoint_labels()

  df <- .make_suppa_multiplicity_magnitude_df(
    suppa_results = suppa_results,
    p_cutoff = p_cutoff
  )

  plot_df <- df |>
    dplyr::filter(n_sig_events <= max_events_to_show)

  global_median_events <- stats::median(plot_df$n_sig_events, na.rm = TRUE)
  global_median_dpsi <- stats::median(plot_df$max_abs_dpsi, na.rm = TRUE)

  global_subtitle <- bquote(
    "Median Events per Gene = " * .(round(global_median_events, 1)) *
      "; Median Max |" * Delta * "PSI| = " * .(round(global_median_dpsi, 2))
  )

  stat_df <- plot_df |>
    dplyr::group_by(timepoint) |>
    dplyr::summarize(
      median_events = stats::median(n_sig_events, na.rm = TRUE),
      mean_events = mean(n_sig_events, na.rm = TRUE),
      median_dpsi = stats::median(max_abs_dpsi, na.rm = TRUE),
      mean_dpsi = mean(max_abs_dpsi, na.rm = TRUE),
      .groups = "drop"
    )

  stat_label_df <- stat_df |>
    dplyr::mutate(
      stat_label = paste0(
        "atop(",
        "'Mean Events = ", round(mean_events, 2), "',",
        "'Mean |' * Delta * 'PSI| = ", round(mean_dpsi, 3), "'",
        ")"
      )
    )

  label_df <- plot_df |>
    dplyr::group_by(timepoint) |>
    dplyr::filter(n_sig_events >= label_min_events) |>
    dplyr::arrange(dplyr::desc(n_sig_events), dplyr::desc(max_abs_dpsi), .by_group = TRUE) |>
    dplyr::slice_head(n = label_top_n_by_events) |>
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
    ggplot2::aes(x = n_sig_events, y = max_abs_dpsi)
  ) +
    ggplot2::geom_vline(
      data = stat_df,
      ggplot2::aes(xintercept = median_events),
      linetype = "dashed",
      linewidth = 0.5,
      color = unname(grey_palette_colors["dark"])
    ) +
    ggplot2::geom_hline(
      data = stat_df,
      ggplot2::aes(yintercept = median_dpsi),
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
      color = splicing_tool_colors[["SUPPA2"]],
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
      labeller = ggplot2::labeller(timepoint = timepoint_labels),
      strip = .suppa_timepoint_strip()
    ) +
    ggplot2::scale_x_continuous(
      limits = c(1, max_events_to_show),
      breaks = 1:max_events_to_show
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2)
    ) +
    ggplot2::labs(
      x = "Significant SUPPA2 Events per Gene",
      y = expression("Max |" * Delta * "PSI| per Gene"),
      title = title,
      subtitle = global_subtitle
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(
        linewidth = 0.25,
        color = unname(grey_palette_colors["light"])
      ),
      panel.grid.major.y = ggplot2::element_line(
        linewidth = 0.25,
        color = unname(grey_palette_colors["light"])
      ),
      axis.text.x = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 12)
    )

  ggplot2::ggsave(
    filename = file.path(outdir, "sup_fig_1_B_suppa_multiplicity_vs_magnitude.pdf"),
    plot = fig_b,
    width = 12.5,
    height = 9
  )

  fig_b
}

# --------------------------------------------------
# Supplementary Figure 1C - SUPPA2 Event Categories
# --------------------------------------------------

plot_suppa_event_categories <- function(
  suppa_results,
  outdir = NULL,
  title = "SUPPA2: Event Categories",
  base_size = 15,
  p_cutoff = 0.10,
  y_max = 0.5
) {
  .check_pkg(c(
    "dplyr",
    "tidyr",
    "ggplot2",
    "scales",
    "ggh4x"
  ))

  outdir <- .set_outdir(outdir, subdir = "sup_fig_1")
  timepoint_order <- .suppa_timepoint_order()
  timepoint_labels <- .suppa_timepoint_labels()

  timepoint_colors <- .suppa_timepoint_colors(shade = "light")

  event_type_labels <- c(
    "AF" = "Alternative\nFirst Exon",
    "MX" = "Mutually\nExclusive Exons",
    "SE" = "Skipped\nExon",
    "AL" = "Alternative\nLast Exon",
    "A5" = "Alternative\n5' Splice Site",
    "A3" = "Alternative\n3' Splice Site"
  )

  event_type_order <- unname(event_type_labels[c("AF", "MX", "SE", "AL", "A5", "A3")])

  category_subset <- category_base_colors[paste("Category", seq_along(event_type_order))]

  category_shades <- .generate_shaded_palette(
    base_colors = category_subset,
    labels = c("light", "base", "dark")
  )

  event_strip_colors <- unname(vapply(
    category_shades,
    function(x) x["light"],
    character(1)
  ))

  df <- .bind_suppa_events_for_plots(
    suppa_results,
    sig_only = TRUE,
    p_cutoff = p_cutoff
  ) |>
    dplyr::filter(
      !is.na(event_type),
      !is.na(pvalue),
      pvalue < p_cutoff,
      event_type %in% names(event_type_labels)
    ) |>
    dplyr::mutate(
      timepoint = factor(as.character(timepoint), levels = timepoint_order),
      event_type_label = dplyr::recode(event_type, !!!event_type_labels),
      event_type_label = factor(event_type_label, levels = event_type_order)
    ) |>
    dplyr::count(timepoint, event_type_label, name = "n_events") |>
    tidyr::complete(
      timepoint = factor(timepoint_order, levels = timepoint_order),
      event_type_label = factor(event_type_order, levels = event_type_order),
      fill = list(n_events = 0)
    ) |>
    dplyr::group_by(timepoint) |>
    dplyr::mutate(
      total_events = sum(n_events),
      prop = dplyr::if_else(total_events > 0, n_events / total_events, 0)
    ) |>
    dplyr::ungroup()

  fig_c <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = timepoint, y = prop, fill = timepoint)
  ) +
    ggplot2::geom_col(width = 0.75) +
    ggh4x::facet_wrap2(
      ~event_type_label,
      nrow = 1,
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = event_strip_colors,
          color = unname(grey_palette_colors["black"])
        ),
        text_x = ggh4x::elem_list_text(
          color = unname(grey_palette_colors["black"]),
          face = "bold",
          size = 12
        )
      )
    ) +
    ggplot2::scale_fill_manual(
      values = timepoint_colors,
      breaks = timepoint_order,
      labels = timepoint_labels
    ) +
    ggplot2::scale_x_discrete(labels = timepoint_labels, drop = FALSE) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::coord_cartesian(ylim = c(0, y_max)) +
    ggplot2::labs(
      x = NULL,
      y = "% Total SUPPA2 Events",
      fill = NULL,
      title = title
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )

  ggplot2::ggsave(
    filename = file.path(outdir, "sup_fig_1_C_suppa_event_categories.pdf"),
    plot = fig_c,
    width = 12,
    height = 9
  )

  fig_c
}

# --------------------------------------------------
# Supplementary Figure 1D - SUPPA2 Differential Splicing Events
# --------------------------------------------------

plot_suppa_volcano <- function(
  suppa_results,
  outdir = NULL,
  title = "SUPPA2: Differential Splicing Events",
  base_size = 15,
  p_cutoff = 0.10,
  dpsi_cutoff = 0.20,
  label_top_n = 15,
  max_neglog10_p = 5
) {
  .check_pkg(c(
    "dplyr",
    "ggplot2",
    "ggrepel",
    "ggh4x",
    "ggrastr"
  ))

  outdir <- .set_outdir(outdir, subdir = "sup_fig_1")

  timepoint_order <- .suppa_timepoint_order()
  timepoint_labels <- .suppa_timepoint_labels()

  sig_class_order <- c(
    "not_sig",
    "high_dpsi",
    "sig",
    "sig_high_dpsi"
  )

  suppa_shades <- .generate_shaded_palette(
    base_colors = splicing_tool_colors["SUPPA2"],
    labels = c("light", "base", "dark")
  )

  volcano_colors <- c(
    not_sig = unname(grey_palette_colors["mid"]),
    high_dpsi = unname(grey_palette_colors["dark"]),
    sig = unname(suppa_shades$SUPPA2["light"]),
    sig_high_dpsi = unname(suppa_shades$SUPPA2["dark"])
  )

  volcano_labels <- as.expression(c(
    "Not Significant",
    bquote("|" * Delta * "PSI| >= " * .(dpsi_cutoff)),
    bquote("FDR < " * .(p_cutoff)),
    bquote("|" * Delta * "PSI| >= " * .(dpsi_cutoff) * " and FDR < " * .(p_cutoff))
  ))

  df <- .bind_suppa_events_for_plots(
    suppa_results,
    sig_only = FALSE,
    p_cutoff = p_cutoff
  ) |>
    dplyr::filter(
      !is.na(dpsi),
      !is.na(pvalue),
      is.finite(dpsi),
      is.finite(pvalue),
      pvalue > 0
    ) |>
    dplyr::mutate(
      timepoint = factor(as.character(timepoint), levels = timepoint_order),
      sig_class = dplyr::case_when(
        pvalue < p_cutoff & abs(dpsi) >= dpsi_cutoff ~ "sig_high_dpsi",
        pvalue < p_cutoff ~ "sig",
        abs(dpsi) >= dpsi_cutoff ~ "high_dpsi",
        TRUE ~ "not_sig"
      ),
      sig_class = factor(sig_class, levels = sig_class_order),
      neglog10_pvalue = pmin(-log10(pvalue), max_neglog10_p)
    )

  label_df <- df |>
    dplyr::filter(
      pvalue < p_cutoff,
      abs(dpsi) >= dpsi_cutoff
    ) |>
    dplyr::group_by(timepoint) |>
    dplyr::arrange(pvalue, dplyr::desc(abs(dpsi)), .by_group = TRUE) |>
    dplyr::slice_head(n = label_top_n) |>
    dplyr::ungroup()

  fig_d <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = dpsi, y = neglog10_pvalue)
  ) +
    ggrastr::geom_point_rast(
      ggplot2::aes(color = sig_class),
      alpha = 0.75,
      size = 1.75,
      raster.dpi = 600
    ) +
    ggplot2::geom_vline(
      xintercept = c(-dpsi_cutoff, dpsi_cutoff),
      linetype = "dashed",
      linewidth = 0.5,
      color = unname(grey_palette_colors["dark"])
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(p_cutoff),
      linetype = "dashed",
      linewidth = 0.5,
      color = unname(grey_palette_colors["dark"])
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(label = symbol),
      size = 3,
      max.overlaps = Inf,
      max.iter = 50000,
      max.time = 5,
      box.padding = 0.75,
      point.padding = 0.75,
      min.segment.length = 0.1,
      segment.alpha = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::coord_flip(
      xlim = c(-1, 1),
      ylim = c(0, max_neglog10_p)
    ) +
    ggh4x::facet_wrap2(
      ~timepoint,
      nrow = 1,
      labeller = ggplot2::labeller(timepoint = timepoint_labels),
      strip = .suppa_timepoint_strip(shade = "light")
    ) +
    ggplot2::scale_color_manual(
      values = volcano_colors,
      breaks = sig_class_order,
      labels = volcano_labels,
      drop = FALSE,
      name = NULL
    ) +
    ggplot2::labs(
      x = expression(Delta * "PSI"),
      y = expression("-log"[10] * "(p-value)"),
      title = title
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )

  ggplot2::ggsave(
    filename = file.path(outdir, "sup_fig_1_D_suppa_volcano_timepoints.pdf"),
    plot = fig_d,
    width = 12,
    height = 9
  )

  fig_d
}

# --------------------------------------------------
# Convenience Wrapper
# --------------------------------------------------

plot_sup_fig_1_all <- function(
  suppa_results,
  outdir = NULL,
  save_tables = FALSE,
  p_cutoff = 0.10,
  dpsi_cutoff = 0.20,
  label_top_n = 15,
  max_neglog10_p = 5,
  max_events_to_show = 10,
  label_top_n_by_events = 5,
  label_min_events = 5
) {
  outdir <- .set_outdir(outdir, subdir = "sup_fig_1")

  sup_fig_1_a <- plot_suppa_overlap(
    suppa_results = suppa_results,
    outdir = outdir
  )

  sup_fig_1_b <- plot_suppa_multiplicity_vs_magnitude(
    suppa_results = suppa_results,
    outdir = outdir,
    p_cutoff = p_cutoff,
    max_events_to_show = max_events_to_show,
    label_top_n_by_events = label_top_n_by_events,
    label_min_events = label_min_events
  )

  sup_fig_1_c <- plot_suppa_event_categories(
    suppa_results = suppa_results,
    outdir = outdir,
    p_cutoff = p_cutoff
  )

  sup_fig_1_d <- plot_suppa_volcano(
    suppa_results = suppa_results,
    outdir = outdir,
    p_cutoff = p_cutoff,
    dpsi_cutoff = dpsi_cutoff,
    label_top_n = label_top_n,
    max_neglog10_p = max_neglog10_p
  )

  invisible(list(
    sup_fig_1_a = sup_fig_1_a,
    sup_fig_1_b = sup_fig_1_b,
    sup_fig_1_c = sup_fig_1_c,
    sup_fig_1_d = sup_fig_1_d
  ))
}
