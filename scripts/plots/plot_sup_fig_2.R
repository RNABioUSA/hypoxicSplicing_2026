# ============================================================
# plot_sup_fig_2.R
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

.dexseq_timepoint_order <- function() {
  c("H1", "H3", "H24")
}

.dexseq_timepoint_labels <- function() {
  stats::setNames(
    names(timepoint_base_colors)[grepl("Hypoxia", names(timepoint_base_colors))],
    .dexseq_timepoint_order()
  )
}

.dexseq_timepoint_colors <- function(shade = "light") {
  timepoints <- .dexseq_timepoint_order()

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

.dexseq_timepoint_strip <- function(shade = "light") {
  .check_pkg("ggh4x")

  ggh4x::strip_themed(
    background_x = ggh4x::elem_list_rect(
      fill = unname(.dexseq_timepoint_colors(shade = shade)),
      color = unname(grey_palette_colors["black"])
    ),
    text_x = ggh4x::elem_list_text(
      color = unname(grey_palette_colors["black"]),
      face = "bold"
    )
  )
}

.dexseq_tool_shades <- function(labels = c("light", "base", "dark")) {
  .generate_shaded_palette(
    base_colors = splicing_tool_colors["DEXSeq"],
    labels = labels
  )$DEXSeq
}

.dexseq_method_colors <- function() {
  cols <- .make_shades_one(
    base_color = splicing_tool_colors[["DEXSeq"]],
    labels = c("summarizeOverlaps", "HTSeq", "featureCounts")
  )
  cols[c("HTSeq", "summarizeOverlaps", "featureCounts")]
}

.dexseq_h3_results_dir <- function() {
  outdir <- file.path(.get_results_dir(), "analysis", "dexseq_h3_figure")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  outdir
}

.bind_dexseq_exons <- function(dexseq_results) {
  .check_pkg("dplyr")

  if (!("results" %in% names(dexseq_results))) {
    stop("Expected dexseq_results$results.")
  }
  if (!("exon_full_by_tp" %in% names(dexseq_results$results))) {
    stop("Expected dexseq_results$results$exon_full_by_tp.")
  }

  timepoint_order <- .dexseq_timepoint_order()

  dplyr::bind_rows(dexseq_results$results$exon_full_by_tp, .id = "tp_name") |>
    dplyr::mutate(
      timepoint = dplyr::if_else(
        is.na(timepoint) | timepoint == "",
        tp_name,
        as.character(timepoint)
      ),
      timepoint = factor(timepoint, levels = timepoint_order),
      symbol = dplyr::if_else(
        is.na(symbol) | symbol == "",
        ensgene,
        symbol
      )
    ) |>
    dplyr::select(-tp_name)
}

.get_dexseq_effect_col <- function(df, tp) {
  candidates <- c(
    paste0("log2fold_", tp, "_C1"),
    "log2fold"
  )

  effect_col <- candidates[candidates %in% names(df)][1]

  if (is.na(effect_col) || length(effect_col) == 0) {
    stop(
      "Could not find DEXSeq effect column for ", tp,
      ". Checked: ", paste(candidates, collapse = ", ")
    )
  }

  effect_col
}

.add_dexseq_feature_id <- function(df) {
  df |>
    dplyr::mutate(
      feature_id_plot = dplyr::case_when(
        "feature_id" %in% names(df) ~ as.character(feature_id),
        "exon_id" %in% names(df) ~ as.character(exon_id),
        TRUE ~ paste0(ensgene, "_exon_", dplyr::row_number())
      )
    )
}

.dexseq_h3_policy_table <- function() {
  x <- expand.grid(
    counting_method = c("summarizeOverlaps", "featureCounts"),
    count_multi_exon_overlaps = c(TRUE, FALSE),
    include_multimappers = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  htseq_row <- data.frame(
    counting_method = "HTSeq",
    count_multi_exon_overlaps = TRUE,
    include_multimappers = FALSE,
    stringsAsFactors = FALSE
  )

  x <- rbind(x, htseq_row)

  x$mode_tag <- vapply(
    seq_len(nrow(x)),
    function(i) {
      .get_dexseq_mode_tag(
        counting_method = x$counting_method[i],
        count_multi_exon_overlaps = x$count_multi_exon_overlaps[i],
        include_multimappers = x$include_multimappers[i]
      )
    },
    character(1)
  )

  x$timepoint <- "H3"
  x$analysis_label <- paste0("H3.", x$mode_tag)
  x$out_rds <- file.path(.dexseq_h3_results_dir(), paste0("dexseq_results.", x$analysis_label, ".rds"))
  x$out_xlsx <- file.path(.dexseq_h3_results_dir(), paste0("dexseq_results.", x$analysis_label, ".xlsx"))

  x$overlap_policy <- ifelse(x$count_multi_exon_overlaps, "multiOverlap", "strictOverlap")
  x$multimap_policy <- ifelse(x$include_multimappers, "all", "unique")
  x$method_label <- x$counting_method
  x$policy_label <- ifelse(
    x$counting_method == "HTSeq",
    "multiOverlap.unique",
    paste0(x$overlap_policy, ".", x$multimap_policy)
  )
  x$plot_label <- ifelse(
    x$counting_method == "HTSeq",
    "HTSeq\nmultiOverlap.unique",
    paste0(x$overlap_policy, "\n", x$multimap_policy)
  )

  x
}

.annotate_dexseq_h3_summary <- function(summary_tbl) {
  if (!("mode_tag" %in% colnames(summary_tbl))) {
    stop("Summary table must contain a 'mode_tag' column.")
  }

  summary_tbl$counting_method <- sub("\\..*", "", summary_tbl$mode_tag)
  summary_tbl$overlap_policy <- NA_character_
  summary_tbl$multimap_policy <- NA_character_

  is_htseq <- summary_tbl$counting_method == "HTSeq"
  summary_tbl$overlap_policy[is_htseq] <- "multiOverlap"
  summary_tbl$multimap_policy[is_htseq] <- "unique"

  is_other <- !is_htseq
  split_tags <- strsplit(summary_tbl$mode_tag[is_other], "\\.")
  summary_tbl$overlap_policy[is_other] <- vapply(split_tags, `[`, character(1), 2)
  summary_tbl$multimap_policy[is_other] <- vapply(split_tags, `[`, character(1), 3)

  numeric_cols <- c("n_tested_genes", "n_sig_genes", "n_tested_exons", "n_sig_exons")
  for (cc in intersect(numeric_cols, colnames(summary_tbl))) {
    summary_tbl[[cc]] <- suppressWarnings(as.numeric(summary_tbl[[cc]]))
  }

  summary_tbl$method_label <- summary_tbl$counting_method
  summary_tbl$policy_label <- paste0(summary_tbl$overlap_policy, ".", summary_tbl$multimap_policy)
  summary_tbl$plot_label <- paste0(summary_tbl$overlap_policy, "\n", summary_tbl$multimap_policy)

  summary_tbl$method_label <- factor(
    summary_tbl$method_label,
    levels = c("HTSeq", "summarizeOverlaps", "featureCounts")
  )

  summary_tbl$counting_method <- factor(
    summary_tbl$counting_method,
    levels = c("HTSeq", "summarizeOverlaps", "featureCounts")
  )

  summary_tbl$plot_label <- factor(
    summary_tbl$plot_label,
    levels = c(
      "multiOverlap\nall",
      "multiOverlap\nunique",
      "strictOverlap\nall",
      "strictOverlap\nunique"
    )
  )

  summary_tbl
}

# --------------------------------------------------
# Analysis / Data Construction
# --------------------------------------------------

.make_dexseq_multiplicity_magnitude_df <- function(
  dexseq_results,
  padj_cutoff = 0.10,
  sig_only = TRUE
) {
  .check_pkg("dplyr")

  timepoint_order <- .dexseq_timepoint_order()
  dex <- .bind_dexseq_exons(dexseq_results)

  if (!("significant" %in% names(dex))) {
    dex <- dex |>
      dplyr::mutate(significant = !is.na(padj) & padj <= padj_cutoff)
  }

  if (isTRUE(sig_only)) {
    dex <- dex |> dplyr::filter(significant %in% TRUE)
  }

  dplyr::bind_rows(lapply(timepoint_order, function(tp) {
    tbl <- dex |> dplyr::filter(as.character(timepoint) == tp)

    if (nrow(tbl) == 0) {
      return(NULL)
    }

    effect_col <- .get_dexseq_effect_col(tbl, tp)

    tbl |>
      dplyr::mutate(
        effect = .data[[effect_col]],
        abs_effect = abs(effect)
      ) |>
      .add_dexseq_feature_id() |>
      dplyr::filter(
        !is.na(ensgene),
        is.finite(abs_effect)
      ) |>
      dplyr::group_by(timepoint, ensgene, symbol) |>
      dplyr::summarize(
        n_sig_exons = dplyr::n_distinct(feature_id_plot),
        max_abs_exon_effect = max(abs_effect, na.rm = TRUE),
        mean_abs_exon_effect = mean(abs_effect, na.rm = TRUE),
        .groups = "drop"
      )
  })) |>
    dplyr::mutate(tool = "DEXSeq")
}

.make_dexseq_direction_consistency_df <- function(
  dexseq_results,
  padj_cutoff = 0.10,
  min_sig_exons = 2
) {
  .check_pkg("dplyr")

  timepoint_order <- .dexseq_timepoint_order()
  dex <- .bind_dexseq_exons(dexseq_results)

  if (!("significant" %in% names(dex))) {
    dex <- dex |>
      dplyr::mutate(significant = !is.na(padj) & padj <= padj_cutoff)
  }

  dplyr::bind_rows(lapply(timepoint_order, function(tp) {
    tbl <- dex |>
      dplyr::filter(as.character(timepoint) == tp, significant %in% TRUE)

    if (nrow(tbl) == 0) {
      return(NULL)
    }

    effect_col <- .get_dexseq_effect_col(tbl, tp)

    tbl |>
      dplyr::mutate(
        effect = .data[[effect_col]],
        effect_direction = dplyr::case_when(
          effect > 0 ~ "positive",
          effect < 0 ~ "negative",
          TRUE ~ NA_character_
        )
      ) |>
      .add_dexseq_feature_id() |>
      dplyr::filter(!is.na(effect_direction), !is.na(ensgene)) |>
      dplyr::group_by(timepoint, ensgene, symbol) |>
      dplyr::summarize(
        n_sig_exons = dplyr::n_distinct(feature_id_plot),
        n_pos = sum(effect_direction == "positive", na.rm = TRUE),
        n_neg = sum(effect_direction == "negative", na.rm = TRUE),
        direction_class = dplyr::case_when(
          n_pos > 0 & n_neg == 0 ~ "All Positive",
          n_neg > 0 & n_pos == 0 ~ "All Negative",
          n_pos > 0 & n_neg > 0 ~ "Mixed Direction",
          TRUE ~ NA_character_
        ),
        .groups = "drop"
      ) |>
      dplyr::filter(n_sig_exons >= min_sig_exons, !is.na(direction_class))
  })) |>
    dplyr::mutate(
      direction_class = factor(
        direction_class,
        levels = c("All Positive", "All Negative", "Mixed Direction")
      )
    )
}

.make_dexseq_direction_enrichment_df <- function(
  dexseq_results,
  padj_cutoff = 0.10,
  min_sig_exons = 2
) {
  .check_pkg(c("dplyr", "tidyr"))

  obs_gene_df <- .make_dexseq_direction_consistency_df(
    dexseq_results = dexseq_results,
    padj_cutoff = padj_cutoff,
    min_sig_exons = min_sig_exons
  )

  obs_df <- obs_gene_df |>
    dplyr::count(timepoint, direction_class, name = "n_obs") |>
    dplyr::group_by(timepoint) |>
    dplyr::mutate(obs_prop = n_obs / sum(n_obs)) |>
    dplyr::ungroup()

  expected_df <- obs_gene_df |>
    dplyr::group_by(timepoint) |>
    dplyr::mutate(
      total_pos = sum(n_pos),
      total_neg = sum(n_neg),
      p_pos = total_pos / (total_pos + total_neg),
      p_neg = total_neg / (total_pos + total_neg)
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      exp_all_positive = p_pos^n_sig_exons,
      exp_all_negative = p_neg^n_sig_exons,
      exp_mixed = 1 - exp_all_positive - exp_all_negative
    ) |>
    dplyr::ungroup() |>
    dplyr::summarize(
      exp_all_positive = mean(exp_all_positive, na.rm = TRUE),
      exp_all_negative = mean(exp_all_negative, na.rm = TRUE),
      exp_mixed = mean(exp_mixed, na.rm = TRUE),
      .by = timepoint
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::starts_with("exp_"),
      names_to = "direction_class",
      values_to = "exp_prop"
    ) |>
    dplyr::mutate(
      direction_class = dplyr::recode(
        direction_class,
        exp_all_positive = "All Positive",
        exp_all_negative = "All Negative",
        exp_mixed = "Mixed Direction"
      )
    )

  obs_df |>
    dplyr::left_join(expected_df, by = c("timepoint", "direction_class")) |>
    dplyr::mutate(
      enrichment = obs_prop / exp_prop,
      log2_enrichment = log2(enrichment),
      direction_class = factor(
        direction_class,
        levels = c("All Positive", "All Negative", "Mixed Direction")
      )
    )
}

.make_dexseq_exon_position_map <- function(dexseq_results) {
  .check_pkg("dplyr")

  dex_all <- .bind_dexseq_exons(dexseq_results)

  dex_all |>
    dplyr::mutate(
      exon_start = as.numeric(genomicData.start),
      exon_end = as.numeric(genomicData.end),
      strand = as.character(genomicData.strand),
      exon_key = paste(ensgene, feature_id, exon_start, exon_end, strand, sep = "|")
    ) |>
    dplyr::distinct(
      ensgene, feature_id, exon_key,
      exon_start, exon_end, strand
    ) |>
    dplyr::filter(
      !is.na(ensgene),
      !is.na(feature_id),
      is.finite(exon_start),
      is.finite(exon_end),
      strand %in% c("+", "-")
    ) |>
    dplyr::mutate(
      exon_order_coord = dplyr::if_else(
        strand == "-",
        -exon_start,
        exon_start
      )
    ) |>
    dplyr::group_by(ensgene) |>
    dplyr::arrange(exon_order_coord, exon_end, .by_group = TRUE) |>
    dplyr::mutate(
      exon_rank = dplyr::row_number(),
      n_exon_bins_gene = dplyr::n(),
      exon_position = dplyr::case_when(
        exon_rank == 1 ~ "First Exon",
        exon_rank == n_exon_bins_gene ~ "Last Exon",
        TRUE ~ "Internal Exon"
      ),
      exon_position = factor(
        exon_position,
        levels = c("First Exon", "Internal Exon", "Last Exon")
      )
    ) |>
    dplyr::ungroup()
}

.make_dexseq_exon_position_enrichment_df <- function(
  dexseq_results,
  padj_cutoff = 0.10
) {
  .check_pkg("dplyr")

  pos_map <- .make_dexseq_exon_position_map(dexseq_results)

  dex_all <- .bind_dexseq_exons(dexseq_results) |>
    dplyr::mutate(
      exon_start = as.numeric(genomicData.start),
      exon_end = as.numeric(genomicData.end),
      strand = as.character(genomicData.strand),
      exon_key = paste(ensgene, feature_id, exon_start, exon_end, strand, sep = "|")
    ) |>
    dplyr::left_join(
      pos_map |> dplyr::select(exon_key, exon_position),
      by = "exon_key"
    ) |>
    dplyr::filter(!is.na(exon_position))

  if (!("significant" %in% names(dex_all))) {
    dex_all <- dex_all |>
      dplyr::mutate(significant = !is.na(padj) & padj <= padj_cutoff)
  }

  bg_df <- dex_all |>
    dplyr::distinct(ensgene, feature_id, exon_start, exon_end, strand, exon_position) |>
    dplyr::count(exon_position, name = "n_bg") |>
    dplyr::mutate(bg_prop = n_bg / sum(n_bg))

  sig_df <- dex_all |>
    dplyr::filter(significant %in% TRUE) |>
    dplyr::count(timepoint, exon_position, name = "n_sig") |>
    dplyr::group_by(timepoint) |>
    dplyr::mutate(sig_prop = n_sig / sum(n_sig)) |>
    dplyr::ungroup()

  sig_df |>
    dplyr::left_join(bg_df, by = "exon_position") |>
    dplyr::mutate(
      enrichment = sig_prop / bg_prop,
      log2_enrichment = log2(enrichment)
    )
}

run_dexseq_h3_strategy_analysis <- function(
  overwrite = FALSE,
  load_existing = FALSE,
  return_results = FALSE
) {
  policy_tbl <- .dexseq_h3_policy_table()

  out <- if (isTRUE(return_results)) {
    vector("list", nrow(policy_tbl))
  } else {
    NULL
  }
  if (isTRUE(return_results)) names(out) <- policy_tbl$mode_tag

  for (i in seq_len(nrow(policy_tbl))) {
    mode_tag <- policy_tbl$mode_tag[i]
    out_rds <- policy_tbl$out_rds[i]
    out_xlsx <- policy_tbl$out_xlsx[i]

    if (file.exists(out_rds) && !overwrite) {
      if (isTRUE(load_existing) && isTRUE(return_results)) {
        message("[DEXSeq Plots] Loading existing: ", mode_tag)
        out[[mode_tag]] <- readRDS(out_rds)
      } else {
        message("[DEXSeq Plots] Skipping existing: ", mode_tag)
      }
      next
    }

    message("[DEXSeq Plots] Running: ", mode_tag)

    res <- run_dexseq(
      timepoints = "H3",
      counting_method = policy_tbl$counting_method[i],
      count_multi_exon_overlaps = policy_tbl$count_multi_exon_overlaps[i],
      include_multimappers = policy_tbl$include_multimappers[i],
      force = overwrite,
      out_rds = out_rds,
      out_xlsx = out_xlsx
    )

    if (isTRUE(return_results)) {
      out[[mode_tag]] <- res
    }

    rm(res)
    gc()
  }

  if (isTRUE(return_results)) {
    attr(out, "policy_table") <- policy_tbl
    return(out)
  }

  invisible(policy_tbl)
}

parse_dexseq_h3_strategy_inputs <- function(
  delete_rds = FALSE,
  overwrite = FALSE
) {
  policy_tbl <- .dexseq_h3_policy_table()
  outdir <- .dexseq_h3_results_dir()

  summary_list <- list()

  for (i in seq_len(nrow(policy_tbl))) {
    mode_tag <- policy_tbl$mode_tag[i]
    rds_file <- policy_tbl$out_rds[i]

    summary_file <- file.path(outdir, paste0("plot_summary.H3.", mode_tag, ".tsv"))
    gene_set_file <- file.path(outdir, paste0("gene_set.H3.", mode_tag, ".txt"))

    if (file.exists(summary_file) && file.exists(gene_set_file) && !overwrite) {
      message("[DEXSeq Plots] Using cached plot inputs: ", mode_tag)
      summary_list[[mode_tag]] <- utils::read.delim(
        summary_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )

      if (isTRUE(delete_rds) && file.exists(rds_file)) {
        file.remove(rds_file)
      }
      next
    }

    if (!file.exists(rds_file)) {
      warning("Missing RDS for mode: ", mode_tag, "\n", rds_file)
      next
    }

    message("[DEXSeq Plots] Extracting plot inputs from RDS: ", mode_tag)
    res <- readRDS(rds_file)

    gene_sig_all <- res$results$gene_sig_all
    exon_sig_all <- res$results$exon_sig_all
    gene_full_all <- res$results$gene_full_all
    exon_full_all <- res$results$exon_full_all

    gene_sig_h3 <- unique(gene_sig_all$gene_id[gene_sig_all$timepoint == "H3"])
    exon_sig_h3 <- unique(exon_sig_all$feature_id[exon_sig_all$timepoint == "H3"])
    gene_full_h3 <- unique(gene_full_all$gene_id[gene_full_all$timepoint == "H3"])
    exon_full_h3 <- unique(exon_full_all$feature_id[exon_full_all$timepoint == "H3"])

    summary_tbl <- data.frame(
      mode_tag = mode_tag,
      counting_method = policy_tbl$counting_method[i],
      overlap_policy = policy_tbl$overlap_policy[i],
      multimap_policy = policy_tbl$multimap_policy[i],
      method_label = policy_tbl$method_label[i],
      policy_label = policy_tbl$policy_label[i],
      n_tested_genes = length(gene_full_h3),
      n_sig_genes = length(gene_sig_h3),
      n_tested_exons = length(exon_full_h3),
      n_sig_exons = length(exon_sig_h3),
      stringsAsFactors = FALSE
    )

    utils::write.table(
      summary_tbl,
      file = summary_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    writeLines(sort(gene_sig_h3), gene_set_file)

    summary_list[[mode_tag]] <- summary_tbl

    rm(res, gene_sig_all, exon_sig_all, gene_full_all, exon_full_all)
    gc()

    if (isTRUE(delete_rds)) {
      file.remove(rds_file)
    }
  }

  if (length(summary_list) == 0) {
    stop("No DEXSeq H3 strategy summaries were available or generated.")
  }

  summary_all <- do.call(rbind, summary_list)
  rownames(summary_all) <- NULL

  combined_summary_file <- file.path(outdir, "plot_summary.H3.all_policies.tsv")
  utils::write.table(
    summary_all,
    file = combined_summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  invisible(summary_all)
}

load_dexseq_h3_strategy_inputs <- function(parse_if_missing = TRUE) {
  policy_tbl <- .dexseq_h3_policy_table()
  outdir <- .dexseq_h3_results_dir()

  summary_file <- file.path(outdir, "plot_summary.H3.all_policies.tsv")
  if (!file.exists(summary_file)) {
    if (isTRUE(parse_if_missing)) {
      summary_tbl <- parse_dexseq_h3_strategy_inputs(delete_rds = FALSE)
    } else {
      stop("Missing summary file: ", summary_file)
    }
  } else {
    summary_tbl <- utils::read.delim(
      summary_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  summary_tbl <- .annotate_dexseq_h3_summary(summary_tbl)

  gene_sets <- lapply(policy_tbl$mode_tag, function(mode_tag) {
    gene_file <- file.path(outdir, paste0("gene_set.H3.", mode_tag, ".txt"))
    if (!file.exists(gene_file)) {
      stop(
        "Missing gene-set file for mode: ", mode_tag,
        "\nRun parse_dexseq_h3_strategy_inputs() first."
      )
    }
    readLines(gene_file, warn = FALSE)
  })
  names(gene_sets) <- policy_tbl$mode_tag

  list(
    summary = summary_tbl,
    gene_sets = gene_sets,
    policy_table = policy_tbl
  )
}

run_dexseq_supplement_analysis <- function(
  dexseq_results,
  outdir = NULL,
  padj_cutoff = 0.10,
  sig_only = TRUE,
  min_sig_exons = 2,
  run_h3_strategy = FALSE,
  parse_h3_strategy = TRUE,
  overwrite_h3_strategy = FALSE,
  delete_h3_rds = FALSE,
  save_tables = FALSE,
  force_recompute = FALSE
) {
  .check_pkg("dplyr")

  outdir <- .set_outdir(outdir, subdir = "sup_fig_2")

  outfile <- file.path(outdir, "sup_fig_2_dexseq_analysis.rds")

  if (!force_recompute && file.exists(outfile)) {
    message("Loading cached DEXSeq analysis from: ", outfile)
    return(readRDS(outfile))
  }

  genes_by_tp <- lapply(.dexseq_timepoint_order(), function(tp) {
    unique(stats::na.omit(dexseq_results$results$gene_sig_by_tp[[tp]]$ensgene))
  })
  names(genes_by_tp) <- .dexseq_timepoint_order()

  multiplicity_df <- .make_dexseq_multiplicity_magnitude_df(
    dexseq_results = dexseq_results,
    padj_cutoff = padj_cutoff,
    sig_only = sig_only
  )

  exon_position_enrichment_df <- .make_dexseq_exon_position_enrichment_df(
    dexseq_results = dexseq_results,
    padj_cutoff = padj_cutoff
  )

  direction_consistency_df <- .make_dexseq_direction_consistency_df(
    dexseq_results = dexseq_results,
    padj_cutoff = padj_cutoff,
    min_sig_exons = min_sig_exons
  )

  direction_enrichment_df <- .make_dexseq_direction_enrichment_df(
    dexseq_results = dexseq_results,
    padj_cutoff = padj_cutoff,
    min_sig_exons = min_sig_exons
  )

  if (isTRUE(run_h3_strategy)) {
    run_dexseq_h3_strategy_analysis(
      overwrite = overwrite_h3_strategy,
      load_existing = FALSE,
      return_results = FALSE
    )
  }

  h3_strategy_inputs <- NULL
  if (isTRUE(parse_h3_strategy)) {
    parse_dexseq_h3_strategy_inputs(
      delete_rds = delete_h3_rds,
      overwrite = overwrite_h3_strategy
    )
    h3_strategy_inputs <- load_dexseq_h3_strategy_inputs(parse_if_missing = FALSE)
  } else {
    h3_strategy_inputs <- tryCatch(
      load_dexseq_h3_strategy_inputs(parse_if_missing = FALSE),
      error = function(e) NULL
    )
  }

  out <- list(
    meta = list(
      outdir = outdir,
      padj_cutoff = padj_cutoff,
      sig_only = sig_only,
      min_sig_exons = min_sig_exons,
      timepoint_order = .dexseq_timepoint_order()
    ),
    overlap = list(
      genes_by_tp = genes_by_tp
    ),
    multiplicity = list(
      gene_df = multiplicity_df
    ),
    exon_direction = list(
      position_enrichment_df = exon_position_enrichment_df,
      direction_consistency_df = direction_consistency_df,
      direction_enrichment_df = direction_enrichment_df
    ),
    h3_strategy = h3_strategy_inputs
  )

  if (isTRUE(save_tables)) {
    utils::write.table(
      multiplicity_df,
      file = file.path(outdir, "sup_fig_2_B_dexseq_multiplicity_magnitude.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    utils::write.table(
      exon_position_enrichment_df,
      file = file.path(outdir, "sup_fig_2_C_dexseq_exon_position_enrichment.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    utils::write.table(
      direction_enrichment_df,
      file = file.path(outdir, "sup_fig_2_C_dexseq_direction_enrichment.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    if (!is.null(h3_strategy_inputs)) {
      utils::write.table(
        h3_strategy_inputs$summary,
        file = file.path(outdir, "sup_fig_2_D_dexseq_h3_strategy_summary.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
    }
  }

  saveRDS(out, file.path(outdir, "sup_fig_2_dexseq_analysis.rds"))
  out
}

# --------------------------------------------------
# Panel A - DEXSeq Gene Overlap Across Timepoints
# --------------------------------------------------

plot_dexseq_overlap_across_timepoints <- function(
  res,
  outdir = NULL,
  title = "DEXSeq: Gene Overlap Across Hypoxia Timepoints"
) {
  .check_pkg("grid")

  outdir <- .set_outdir(outdir, subdir = "sup_fig_2")

  if (is.null(res$overlap$genes_by_tp)) {
    stop("Expected res$overlap$genes_by_tp.")
  }

  fig_a <- plot_feature_overlap_across_timepoints(
    features_by_tp = res$overlap$genes_by_tp,
    title = title,
    shade_label = "TX"
  )

  grid::grid.newpage()
  grid::grid.draw(fig_a$plot)

  grDevices::pdf(
    file.path(outdir, "sup_fig_2_A_dexseq_gene_overlap.pdf"),
    width = 6.5,
    height = 6.5
  )
  grid::grid.newpage()
  grid::grid.draw(fig_a$plot)
  invisible(grDevices::dev.off())

  invisible(fig_a)
}

# --------------------------------------------------
# Panel B - DEXSeq Multiplicity and Magnitude
# --------------------------------------------------

plot_dexseq_multiplicity_vs_magnitude <- function(
  res,
  outdir = NULL,
  title = "DEXSeq: Multiplicity and Magnitude",
  base_size = 15,
  max_features_to_show = 10,
  max_effect_to_show = 5,
  label_top_n_by_features = 5,
  label_min_features = 5,
  point_alpha = 0.75,
  point_size = 3
) {
  .check_pkg(c("dplyr", "ggplot2", "ggrepel", "ggh4x", "scales"))

  outdir <- .set_outdir(outdir, subdir = "sup_fig_2")

  df <- res$multiplicity$gene_df

  plot_df <- df |>
    dplyr::filter(
      n_sig_exons <= max_features_to_show,
      max_abs_exon_effect <= max_effect_to_show
    )

  excluded_df <- df |>
    dplyr::filter(
      n_sig_exons > max_features_to_show |
        max_abs_exon_effect > max_effect_to_show
    )

  global_median_features <- stats::median(plot_df$n_sig_exons, na.rm = TRUE)
  global_median_effect <- stats::median(plot_df$max_abs_exon_effect, na.rm = TRUE)

  global_subtitle <- bquote(
    "Median Exons per Gene = " * .(round(global_median_features, 1)) *
      "; Median Max |Exon log"[2] * "FC| = " * .(round(global_median_effect, 2))
  )

  stat_df <- plot_df |>
    dplyr::group_by(timepoint) |>
    dplyr::summarize(
      median_features = stats::median(n_sig_exons, na.rm = TRUE),
      mean_features = mean(n_sig_exons, na.rm = TRUE),
      median_effect = stats::median(max_abs_exon_effect, na.rm = TRUE),
      mean_effect = mean(max_abs_exon_effect, na.rm = TRUE),
      n_genes = dplyr::n(),
      .groups = "drop"
    )

  stat_label_df <- stat_df |>
    dplyr::mutate(
      stat_label = paste0(
        "Mean Exons = ", round(mean_features, 2),
        "\nMean |Exon log2FC| = ", round(mean_effect, 3)
      )
    )

  label_df <- plot_df |>
    dplyr::group_by(timepoint) |>
    dplyr::filter(n_sig_exons >= label_min_features) |>
    dplyr::arrange(dplyr::desc(n_sig_exons), dplyr::desc(max_abs_exon_effect), .by_group = TRUE) |>
    dplyr::slice_head(n = label_top_n_by_features) |>
    dplyr::ungroup() |>
    dplyr::mutate(label_symbol = symbol)

  repel_df <- plot_df |>
    dplyr::left_join(
      label_df |> dplyr::select(timepoint, ensgene, label_symbol),
      by = c("timepoint", "ensgene")
    )

  fig_b <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = n_sig_exons, y = max_abs_exon_effect)
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
      alpha = 0.9
    ) +
    ggplot2::geom_point(
      color = splicing_tool_colors[["DEXSeq"]],
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
      labeller = ggplot2::labeller(timepoint = .dexseq_timepoint_labels()),
      strip = .dexseq_timepoint_strip(shade = "light")
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
      x = "Significant DEXSeq Exons per Gene",
      y = expression("Max |Exon log"[2] * "FC| per Gene"),
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
    filename = file.path(outdir, "sup_fig_2_B_dexseq_multiplicity_vs_magnitude.pdf"),
    plot = fig_b,
    width = 12.5,
    height = 9
  )

  fig_b
}

# --------------------------------------------------
# Panel C - Counting Method / Strategy
# --------------------------------------------------

.plot_dexseq_h3_strategy_counts <- function(
  res,
  title = "DEXSeq: Counting Method / Strategy Counts",
  base_size = 15
) {
  .check_pkg(c("ggplot2", "ggh4x", "ggpattern"))

  if (is.null(res$h3_strategy)) {
    stop("Expected res$h3_strategy. Run/load H3 strategy inputs first.")
  }

  df <- res$h3_strategy$summary
  dexseq_method_cols <- .dexseq_method_colors()

  df$counting_method <- factor(
    df$counting_method,
    levels = c("HTSeq", "summarizeOverlaps", "featureCounts")
  )

  df$overlap_label <- ifelse(
    df$overlap_policy == "multiOverlap",
    "Multiple Exon\nOverlap",
    "Single Exon\nOverlap"
  )

  df$overlap_label <- factor(
    df$overlap_label,
    levels = c("Multiple Exon\nOverlap", "Single Exon\nOverlap")
  )

  df$alignment_label <- ifelse(
    df$multimap_policy == "all",
    "All Alignments",
    "Unique Alignments"
  )

  df$alignment_label <- factor(
    df$alignment_label,
    levels = c("Unique Alignments", "All Alignments")
  )

  df <- df[
    !is.na(df$counting_method) &
      !is.na(df$overlap_label) &
      !is.na(df$alignment_label) &
      !is.na(df$n_sig_genes), ,
    drop = FALSE
  ]

  stripe_col_light <- .make_shades_one(
    base_color = dexseq_method_cols[["summarizeOverlaps"]],
    labels = c("light_light", "light_base", "light_dark")
  )[["light_light"]]

  stripe_col_dark <- .make_shades_one(
    base_color = dexseq_method_cols[["featureCounts"]],
    labels = c("dark_light", "dark_base", "dark_dark")
  )[["dark_dark"]]

  strip_text_cols <- c(
    "HTSeq" = stripe_col_light,
    "summarizeOverlaps" = stripe_col_dark,
    "featureCounts" = stripe_col_light
  )

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = overlap_label,
      y = n_sig_genes,
      fill = counting_method,
      pattern = alignment_label
    )
  ) +
    ggpattern::geom_col_pattern(
      position = ggplot2::position_dodge(width = 0.75),
      width = 0.55,
      color = unname(grey_palette_colors["black"]),
      pattern_colour = stripe_col_dark,
      pattern_fill = stripe_col_dark,
      pattern_angle = 45,
      pattern_density = 0.18,
      pattern_spacing = 0.035,
      pattern_key_scale_factor = 0.7
    ) +
    ggh4x::facet_wrap2(
      ~counting_method,
      nrow = 1,
      scales = "free_x",
      strip = ggh4x::strip_themed(
        background_x = ggh4x::elem_list_rect(
          fill = dexseq_method_cols[c("HTSeq", "summarizeOverlaps", "featureCounts")]
        ),
        text_x = ggh4x::elem_list_text(
          colour = strip_text_cols[c("HTSeq", "summarizeOverlaps", "featureCounts")],
          face = "bold"
        )
      )
    ) +
    ggplot2::scale_fill_manual(
      values = dexseq_method_cols,
      guide = "none"
    ) +
    ggpattern::scale_pattern_manual(
      values = c(
        "All Alignments" = "stripe",
        "Unique Alignments" = "none"
      ),
      name = "Alignment Strategy"
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Significant Genes",
      title = title
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 0.5
      ),
      panel.grid.major = ggplot2::element_line(color = unname(grey_palette_colors["light"])),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 18, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 15),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "right"
    ) +
    ggplot2::guides(
      pattern = ggplot2::guide_legend(
        title = "Alignment Strategy",
        override.aes = list(
          fill = "white",
          colour = "black",
          pattern_colour = stripe_col_dark,
          pattern_fill = stripe_col_dark
        )
      )
    )
}

.plot_dexseq_h3_strategy_upset <- function(
  res,
  title = "DEXSeq: Counting Method / Strategy Overlap",
  base_size = 15,
  min_size = 25
) {
  .check_pkg(c("ComplexUpset", "ggplot2"))

  if (is.null(res$h3_strategy)) {
    stop("Expected res$h3_strategy. Run/load H3 strategy inputs first.")
  }

  dexseq_method_cols <- .dexseq_method_colors()

  bar_col <- .make_shades_one(
    base_color = dexseq_method_cols[["featureCounts"]],
    labels = c("dark_light", "dark_base", "dark_dark")
  )[["dark_dark"]]

  inactive_dot_col <- .make_shades_one(
    base_color = dexseq_method_cols[["summarizeOverlaps"]],
    labels = c("light_light", "light_base", "light_dark")
  )[["light_light"]]

  upset_tags <- c(
    "HTSeq.multiOverlap.unique",
    "summarizeOverlaps.multiOverlap.all",
    "summarizeOverlaps.multiOverlap.unique",
    "summarizeOverlaps.strictOverlap.all",
    "summarizeOverlaps.strictOverlap.unique",
    "featureCounts.multiOverlap.all",
    "featureCounts.multiOverlap.unique",
    "featureCounts.strictOverlap.all",
    "featureCounts.strictOverlap.unique"
  )

  pretty_names <- rev(c(
    "FC: Single / Unique",
    "FC: Single / All",
    "FC: Multi / Unique",
    "FC: Multi / All",
    "SO: Single / Unique",
    "SO: Single / All",
    "SO: Multi / Unique",
    "SO: Multi / All",
    "HT: Multi / Unique"
  ))

  set_metadata <- data.frame(
    set = pretty_names,
    method = c(
      "HTSeq",
      "summarizeOverlaps", "summarizeOverlaps", "summarizeOverlaps", "summarizeOverlaps",
      "featureCounts", "featureCounts", "featureCounts", "featureCounts"
    ),
    stringsAsFactors = FALSE
  )

  set_metadata$method <- factor(
    set_metadata$method,
    levels = c("HTSeq", "summarizeOverlaps", "featureCounts")
  )

  gene_sets <- res$h3_strategy$gene_sets[upset_tags]
  names(gene_sets) <- pretty_names

  all_genes <- sort(unique(unlist(gene_sets, use.names = FALSE)))

  mat <- data.frame(gene_id = all_genes, stringsAsFactors = FALSE)
  for (set_name in pretty_names) {
    mat[[set_name]] <- all_genes %in% gene_sets[[set_name]]
  }

  suppressWarnings(ComplexUpset::upset(
    mat,
    intersect = rev(pretty_names),
    name = "Counting Strategy Overlap",
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
            values = c("bars_color" = bar_col)
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
      mapping = ggplot2::aes(color = method),
      data = set_metadata,
      colors = dexseq_method_cols
    ),
    themes = ComplexUpset::upset_modify_themes(
      list(
        "Significant Genes" = ggplot2::theme(
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_text(size = 15),
          panel.grid.major = ggplot2::element_line(color = unname(grey_palette_colors["light"])),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_rect(fill = "white", colour = NA),
          plot.background = ggplot2::element_rect(fill = "white", colour = NA)
        )
      )
    )
  ) +
    ggplot2::labs(
      title = title,
      color = "Counting Method"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 18, face = "bold"),
      text = ggplot2::element_text(size = base_size)
    ))
}

plot_dexseq_counting_strategy <- function(
  res,
  outdir = NULL,
  base_size = 15,
  min_size = 25,
  rel_heights = c(1, 2)
) {
  .check_pkg(c("patchwork", "ggplot2"))

  outdir <- .set_outdir(outdir, subdir = "sup_fig_2")

  p_counts <- .plot_dexseq_h3_strategy_counts(
    res = res,
    base_size = base_size
  )

  p_upset <- .plot_dexseq_h3_strategy_upset(
    res = res,
    base_size = base_size,
    min_size = min_size
  )

  fig_d <- p_counts / p_upset +
    patchwork::plot_layout(heights = rel_heights)

  ggplot2::ggsave(
    filename = file.path(outdir, "sup_fig_2_C_dexseq_counting_strategy.pdf"),
    plot = fig_d,
    width = 12,
    height = 12
  )

  invisible(list(
    plot = fig_d,
    counts = p_counts,
    upset = p_upset
  ))
}

# --------------------------------------------------
# Panel D - Exon Position / Direction Bias
# --------------------------------------------------

.plot_dexseq_exon_position_enrichment <- function(
  res,
  title = "DEXSeq: Exon Position Enrichment",
  base_size = 15,
  use_log2 = TRUE
) {
  .check_pkg(c("ggplot2", "scales"))

  df <- res$exon_direction$position_enrichment_df

  y_var <- if (isTRUE(use_log2)) "log2_enrichment" else "enrichment"
  y_lab <- if (isTRUE(use_log2)) {
    "Enrichment (log2)"
  } else {
    "Enrichment"
  }
  ref_line <- if (isTRUE(use_log2)) 0 else 1

  position_levels <- c("First Exon", "Internal Exon", "Last Exon")

  position_shades <- .generate_shaded_palette(
    base_colors = category_base_colors[paste("Category", 1:3)],
    labels = c("light", "base", "dark")
  )

  position_colors <- stats::setNames(
    vapply(position_shades, function(x) unname(x["light"]), character(1)),
    position_levels
  )

  df$exon_position <- factor(df$exon_position, levels = position_levels)

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = timepoint, y = .data[[y_var]], fill = exon_position)
  ) +
    ggplot2::geom_hline(
      yintercept = ref_line,
      linetype = "dashed",
      linewidth = 0.4,
      color = unname(grey_palette_colors["dark"])
    ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.75),
      width = 0.65
    ) +
    ggplot2::scale_fill_manual(
      values = position_colors,
      breaks = position_levels,
      drop = FALSE
    ) +
    ggplot2::scale_x_discrete(labels = .dexseq_timepoint_labels(), drop = FALSE) +
    ggplot2::labs(
      x = NULL,
      y = y_lab,
      fill = NULL,
      title = title,
      subtitle = "Observed exon positions relative to the distribution of all tested exon bins"
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10)
    )
}

.plot_dexseq_direction_enrichment <- function(
  res,
  title = "DEXSeq: Exon Direction Pattern Enrichment",
  base_size = 15,
  use_log2 = TRUE,
  min_sig_exons = NULL
) {
  .check_pkg(c("ggplot2", "scales"))

  df <- res$exon_direction$direction_enrichment_df
  if (is.null(min_sig_exons)) min_sig_exons <- res$meta$min_sig_exons

  y_var <- if (isTRUE(use_log2)) "log2_enrichment" else "enrichment"
  y_lab <- if (isTRUE(use_log2)) {
    "Enrichment (log2)"
  } else {
    "Enrichment"
  }
  ref_line <- if (isTRUE(use_log2)) 0 else 1

  direction_levels <- c("All Positive", "All Negative", "Mixed Direction")

  direction_shades <- .generate_shaded_palette(
    base_colors = category_base_colors[paste("Category", 4:6)],
    labels = c("light", "base", "dark")
  )

  direction_colors <- stats::setNames(
    vapply(direction_shades, function(x) unname(x["light"]), character(1)),
    direction_levels
  )

  df$direction_class <- factor(df$direction_class, levels = direction_levels)

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = timepoint, y = .data[[y_var]], fill = direction_class)
  ) +
    ggplot2::geom_hline(
      yintercept = ref_line,
      linetype = "dashed",
      linewidth = 0.4,
      color = unname(grey_palette_colors["dark"])
    ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.75),
      width = 0.65
    ) +
    ggplot2::scale_fill_manual(
      values = direction_colors,
      breaks = direction_levels,
      drop = FALSE
    ) +
    ggplot2::scale_x_discrete(labels = .dexseq_timepoint_labels(), drop = FALSE) +
    ggplot2::labs(
      x = NULL,
      y = y_lab,
      fill = NULL,
      title = title,
      subtitle = paste0(
        "Observed exon direction patterns relative to expectation based on global exon-level direction proportions for genes with >= ",
        min_sig_exons,
        " significant exon bins"
      )
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10)
    )
}

plot_dexseq_exon_direction_bias <- function(
  res,
  outdir = NULL,
  base_size = 15,
  use_log2 = TRUE,
  rel_heights = c(1, 1)
) {
  .check_pkg(c("patchwork", "ggplot2"))

  outdir <- .set_outdir(outdir, subdir = "sup_fig_2")

  p_position <- .plot_dexseq_exon_position_enrichment(
    res = res,
    base_size = base_size,
    use_log2 = use_log2
  )

  p_direction <- .plot_dexseq_direction_enrichment(
    res = res,
    base_size = base_size,
    use_log2 = use_log2
  )

  fig_c <- p_position / p_direction +
    patchwork::plot_layout(
      heights = rel_heights,
      guides = "keep"
    ) &
    ggplot2::theme(legend.position = "bottom")

  ggplot2::ggsave(
    filename = file.path(outdir, "sup_fig_2_D_dexseq_exon_direction_bias.pdf"),
    plot = fig_c,
    width = 12,
    height = 12
  )

  invisible(list(
    plot = fig_c,
    position = p_position,
    direction = p_direction
  ))
}

# --------------------------------------------------
# Convenience Wrapper
# --------------------------------------------------

plot_sup_fig_2_all <- function(
  dexseq_results,
  outdir = NULL,
  padj_cutoff = 0.10,
  sig_only = TRUE,
  min_sig_exons = 2,
  run_h3_strategy = FALSE,
  parse_h3_strategy = TRUE,
  overwrite_h3_strategy = FALSE,
  delete_h3_rds = FALSE,
  save_tables = FALSE,
  max_features_to_show = 10,
  max_effect_to_show = 5,
  label_top_n_by_features = 5,
  label_min_features = 5,
  use_log2 = TRUE,
  strategy_min_size = 25
) {
  outdir <- .set_outdir(outdir, subdir = "sup_fig_2")

  res <- run_dexseq_supplement_analysis(
    dexseq_results = dexseq_results,
    outdir = outdir,
    padj_cutoff = padj_cutoff,
    sig_only = sig_only,
    min_sig_exons = min_sig_exons,
    run_h3_strategy = run_h3_strategy,
    parse_h3_strategy = parse_h3_strategy,
    overwrite_h3_strategy = overwrite_h3_strategy,
    delete_h3_rds = delete_h3_rds,
    save_tables = save_tables
  )

  fig_a <- plot_dexseq_overlap_across_timepoints(
    res = res,
    outdir = outdir
  )

  fig_b <- plot_dexseq_multiplicity_vs_magnitude(
    res = res,
    outdir = outdir,
    max_features_to_show = max_features_to_show,
    max_effect_to_show = max_effect_to_show,
    label_top_n_by_features = label_top_n_by_features,
    label_min_features = label_min_features
  )

  fig_c <- plot_dexseq_exon_direction_bias(
    res = res,
    outdir = outdir,
    use_log2 = use_log2
  )

  fig_d <- NULL
  if (!is.null(res$h3_strategy)) {
    fig_d <- plot_dexseq_counting_strategy(
      res = res,
      outdir = outdir,
      min_size = strategy_min_size
    )
  } else {
    warning("Skipping Panel D because H3 strategy inputs are unavailable.")
  }

  invisible(list(
    analysis = res,
    fig3_a = fig_a,
    fig3_b = fig_b,
    fig3_c = fig_c,
    fig3_d = fig_d
  ))
}
