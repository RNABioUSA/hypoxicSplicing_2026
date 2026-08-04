# scripts/utils/color_palette.R

# ============================================
# Base Colors
# ============================================

timepoint_base_colors <- c(
  "Normoxia"      = "#C75B5C",
  "Hypoxia (1H)"  = "#5F87A8",
  "Hypoxia (3H)"  = "#7F9F55",
  "Hypoxia (24H)" = "#8A7FAE"
)

splicing_tool_colors <- c(
  "DEXSeq"  = "#C9877A",
  "DRIMSeq" = "#E2C77C",
  "SUPPA2"  = "#7FA7B2",
  "Other"   = "#8FA86A"
)

category_base_colors <- c(
  "Category 1" = "#C57A74",
  "Category 2" = "#CE865F",
  "Category 3" = "#D8AE63",
  "Category 4" = "#4D8A88",
  "Category 5" = "#9B5F6F",
  "Category 6" = "#8A735C"
)


grey_palette_colors <- c(
  light = "#E8E8E8",
  mid   = "#B5B5B5",
  dark  = "#7A7A7A",
  black = "#4A4A4A"
)

# ============================================
# Heatmap Colors
# ============================================

heat_colors <- c(
  "#C45A52",
  "#D98A66",
  "#E6B97A",
  "#F0D8A5",
  "#F3F0D4",
  "#D7E6EB",
  "#A9C7D4",
  "#7FA5C0",
  "#587EA2"
)

# ============================================
# Dot Plot Colors
# ============================================

fgsea_colors_up <- c(
  "#F0D8D8",
  "#E2AAAA",
  "#C75B5C",
  "#974546"
)

fgsea_colors_down <- c(
  "#D7E6EB",
  "#A9C7D4",
  "#5F87A8",
  "#486780"
)

# ============================================
# Generic Shade Helpers
# ============================================

.blend_color <- function(col, target = c("white", "black"), amount) {
  target <- match.arg(target)

  rgb_mat <- grDevices::col2rgb(col) / 255

  target_rgb <- switch(target,
    white = matrix(c(1, 1, 1), nrow = 3),
    black = matrix(c(0, 0, 0), nrow = 3)
  )

  out <- rgb_mat * (1 - amount) + target_rgb * amount
  grDevices::rgb(out[1], out[2], out[3])
}

.make_shades_one <- function(base_color,
                             labels,
                             light_amount = 0.48,
                             dark_amount = 0.24) {
  if (length(labels) != 3) {
    stop("labels must have length 3.")
  }

  stats::setNames(
    c(
      .blend_color(base_color, "white", light_amount),
      base_color,
      .blend_color(base_color, "black", dark_amount)
    ),
    labels
  )
}

.generate_shaded_palette <- function(base_colors,
                                     labels = c("light", "base", "dark"),
                                     light_amount = 0.48,
                                     dark_amount = 0.24) {
  if (is.null(names(base_colors)) || any(names(base_colors) == "")) {
    stop("base_colors must be a named vector.")
  }

  out <- lapply(
    base_colors,
    .make_shades_one,
    labels = labels,
    light_amount = light_amount,
    dark_amount = dark_amount
  )

  names(out) <- names(base_colors)
  out
}

# ============================================
# Timepoint Level Shade Helpers
# ============================================

.get_condition_colors <- function(conditions) {
  color_map <- c(
    "C1"  = unname(timepoint_base_colors["Normoxia"]),
    "H1"  = unname(timepoint_base_colors["Hypoxia (1H)"]),
    "H3"  = unname(timepoint_base_colors["Hypoxia (3H)"]),
    "H24" = unname(timepoint_base_colors["Hypoxia (24H)"])
  )

  cols <- color_map[conditions]

  if (any(is.na(cols))) {
    stop("Missing condition colors for: ", paste(conditions[is.na(cols)], collapse = ", "))
  }
  cols
}

.get_condition_shaded_colors <- function(
  conditions,
  shade = "light",
  light_amount = 0.48,
  dark_amount = 0.24
) {
  color_map <- c(
    "C1"  = unname(timepoint_base_colors["Normoxia"]),
    "H1"  = unname(timepoint_base_colors["Hypoxia (1H)"]),
    "H3"  = unname(timepoint_base_colors["Hypoxia (3H)"]),
    "H24" = unname(timepoint_base_colors["Hypoxia (24H)"])
  )

  base_cols <- color_map[conditions]

  if (any(is.na(base_cols))) {
    stop("Missing condition colors for: ", paste(conditions[is.na(base_cols)], collapse = ", "))
  }

  shaded <- .generate_shaded_palette(
    base_colors = base_cols,
    labels = c("light", "base", "dark"),
    light_amount = light_amount,
    dark_amount = dark_amount
  )

  stats::setNames(
    vapply(shaded, function(x) unname(x[[shade]]), character(1)),
    names(shaded)
  )
}

.get_timepoint_shaded_colors <- function(
  timepoints = c("H1", "H3", "H24"),
  labels = c("TX", "OVERLAP", "GENE")
) {
  tp_label_map <- c(
    H1 = "Hypoxia (1H)",
    H3 = "Hypoxia (3H)",
    H24 = "Hypoxia (24H)"
  )

  missing_tps <- setdiff(timepoints, names(tp_label_map))
  if (length(missing_tps) > 0) {
    stop("Unsupported timepoints: ", paste(missing_tps, collapse = ", "))
  }

  tp_base_colors <- stats::setNames(
    unname(timepoint_base_colors[tp_label_map[timepoints]]),
    timepoints
  )

  .generate_shaded_palette(
    base_colors = tp_base_colors,
    labels = labels
  )
}

.get_timepoint_level_colors <- function(
  timepoints = c("H1", "H3", "H24"),
  light_label = "TX",
  base_label = "OVERLAP",
  dark_label = "GENE"
) {
  shaded <- .get_timepoint_shaded_colors(
    timepoints = timepoints,
    labels = c(light_label, base_label, dark_label)
  )

  out <- c()
  for (tp in timepoints) {
    out[paste0(tp, "_", light_label)] <- shaded[[tp]][[light_label]]
    out[paste0(tp, "_", base_label)] <- shaded[[tp]][[base_label]]
    out[paste0(tp, "_", dark_label)] <- shaded[[tp]][[dark_label]]
  }

  out
}
