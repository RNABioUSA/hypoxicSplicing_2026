#!/usr/bin/env Rscript

# ==============================================================================
# make_supplementary_tables_v2.R
# ------------------------------------------------------------------------------
# Builds a manuscript-ready Excel workbook of curated supplementary tables for
# the hypoxic splicing manuscript.
#
# Expected project layout:
#   $PROJECT_ROOT/results/analysis/deseq/*.rds
#   $PROJECT_ROOT/results/analysis/dexseq/*.rds
#   $PROJECT_ROOT/results/analysis/drimseq/*.rds
#   $PROJECT_ROOT/results/analysis/suppa/standard/C1/*.rds
#   $PROJECT_ROOT/results/plots/fig_1/*.rds
#   $PROJECT_ROOT/results/plots/fig_4/*.rds
#   $PROJECT_ROOT/results/plots/fig_5/*.rds
#   $PROJECT_ROOT/results/plots/fig_6/*.rds
#
# Usage:
#   export PROJECT_ROOT=/path/to/repo
#   Rscript scripts/make_supplementary_tables.R
#
# Optional:
#   Rscript scripts/make_supplementary_tables.R /path/to/output.xlsx
#
# Notes:
#   - This workbook contains significant method-level results and the COMPLETE
#     significant Figure 1 cluster-level GO enrichment output (not only the
#     15 representative terms displayed in Figure 1C).
#   - Full unfiltered RDS objects remain available in the analysis output folders.
#   - DRIMSeq and IsoformSwitchAnalyzeR transcript results are annotated from the
#     same project GENCODE GTF used for quantification (GENCODE v45).
# ==============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required. Install it with install.packages('openxlsx').")
  }
})

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
if (!nzchar(PROJECT_ROOT)) {
  stop("PROJECT_ROOT environment variable is not set. Example: export PROJECT_ROOT=/path/to/repo")
}

args <- commandArgs(trailingOnly = TRUE)
OUT_XLSX <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path(PROJECT_ROOT, "results", "supplementary_data", "Supplementary_Data.xlsx")
}

timepoints <- c("H1", "H3", "H24")
timepoint_labels <- c(H1 = "1 hour hypoxia", H3 = "3 hours hypoxia", H24 = "24 hours hypoxia")

# ------------------------------------------------------------------------------
# Optional project helpers
# ------------------------------------------------------------------------------
HELPERS_FILE <- file.path(PROJECT_ROOT, "scripts", "utils", "helpers.R")
if (file.exists(HELPERS_FILE)) {
  source(HELPERS_FILE)
}

strip_ens_version <- function(x) {
  if (exists(".strip_ens_version", mode = "function")) {
    return(.strip_ens_version(x))
  }
  sub("\\.[0-9]+$", "", as.character(x))
}

# ------------------------------------------------------------------------------
# General helpers
# ------------------------------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x

message_step <- function(...) message("[supplement] ", ...)

first_existing <- function(paths, required = FALSE, label = "file") {
  hit <- paths[file.exists(paths)]
  if (length(hit) > 0) {
    return(hit[[1]])
  }
  if (isTRUE(required)) {
    stop("Could not locate ", label, ". Tried:\n  - ", paste(paths, collapse = "\n  - "))
  }
  NA_character_
}

# Reuse the Figure 6 extraction and GENCODE annotation functions so the figure
# and supplemental workbook cannot silently diverge.
FIG6_SCRIPT <- first_existing(
  c(
    Sys.getenv("PLOT_FIG6_SCRIPT"),
    file.path(PROJECT_ROOT, "scripts", "plots", "plot_fig_6.R"),
    file.path(PROJECT_ROOT, "scripts", "plot_fig_6.R")
  ),
  required = TRUE,
  label = "plot_fig_6.R"
)
source(FIG6_SCRIPT)

load_rds_optional <- function(paths, label) {
  path <- first_existing(paths, required = FALSE, label = label)
  if (is.na(path)) {
    warning("Skipping ", label, "; no RDS found.")
    return(NULL)
  }
  message_step("Loading ", label, ": ", path)
  readRDS(path)
}

clean_for_excel <- function(df) {
  if (is.null(df)) {
    return(data.frame())
  }
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(df) == 0 && ncol(df) == 0) {
    return(data.frame(note = "No records available."))
  }

  # Convert list columns to compact strings.
  for (nm in names(df)) {
    if (is.list(df[[nm]]) && !is.data.frame(df[[nm]])) {
      df[[nm]] <- vapply(df[[nm]], function(z) paste(z, collapse = ";"), character(1))
    }
  }

  # Avoid Excel treating long IDs as formulas or scientific notation where possible.
  char_cols <- vapply(df, is.character, logical(1))
  df[char_cols] <- lapply(df[char_cols], function(x) {
    x <- ifelse(is.na(x), NA_character_, x)
    x <- gsub("[\r\n]+", " ", x)
    x
  })

  df
}

safe_sheet_name <- function(x) {
  x <- gsub("[][\\:*?/\\]", "_", x)
  substr(x, 1, 31)
}

add_sheet <- function(wb, sheet_name, df, freeze_col = 1, freeze_row = 1) {
  sheet_name <- safe_sheet_name(sheet_name)
  if (sheet_name %in% names(wb)) openxlsx::removeWorksheet(wb, sheet_name)
  openxlsx::addWorksheet(wb, sheet_name, gridLines = FALSE)

  df <- clean_for_excel(df)
  openxlsx::writeData(wb, sheet_name, df, withFilter = nrow(df) > 0)

  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    fgFill = "#1F4E78",
    fontColour = "#FFFFFF",
    halign = "center",
    valign = "center",
    border = "Bottom"
  )
  body_style <- openxlsx::createStyle(valign = "top")
  wrap_style <- openxlsx::createStyle(wrapText = TRUE, valign = "top")

  if (ncol(df) > 0) {
    openxlsx::addStyle(wb, sheet_name, header_style, rows = 1, cols = seq_len(ncol(df)), gridExpand = TRUE)
    if (nrow(df) > 0) {
      openxlsx::addStyle(wb, sheet_name, body_style, rows = 2:(nrow(df) + 1), cols = seq_len(ncol(df)), gridExpand = TRUE, stack = TRUE)
    }

    # Widths: cap text-heavy columns to keep workbook readable.
    widths <- pmin(pmax(nchar(names(df)) + 2, 10), 28)
    text_heavy <- grepl("gene|genes|description|term|label|note|README|pathway", names(df), ignore.case = TRUE)
    widths[text_heavy] <- pmin(pmax(widths[text_heavy], 18), 45)
    openxlsx::setColWidths(wb, sheet_name, cols = seq_len(ncol(df)), widths = widths)
    if (any(text_heavy)) {
      openxlsx::addStyle(wb, sheet_name, wrap_style, rows = 1:(nrow(df) + 1), cols = which(text_heavy), gridExpand = TRUE, stack = TRUE)
    }
  }

  openxlsx::freezePane(wb, sheet_name, firstActiveRow = freeze_row + 1, firstActiveCol = freeze_col + 1)
  invisible(sheet_name)
}

select_existing <- function(df, cols) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  cols <- cols[cols %in% colnames(df)]
  if (length(cols) == 0) {
    return(df)
  }
  df[, unique(cols), drop = FALSE]
}

add_missing_columns <- function(df, cols) {
  for (nm in cols) if (!nm %in% colnames(df)) df[[nm]] <- NA
  df
}

# rbind helper that preserves all columns across data frames.
rbind_fill_base <- function(dfs) {
  dfs <- Filter(function(x) !is.null(x) && is.data.frame(x), dfs)
  if (length(dfs) == 0) {
    return(data.frame())
  }
  cols <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  dfs <- lapply(dfs, function(x) {
    missing <- setdiff(cols, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[, cols, drop = FALSE]
  })
  do.call(rbind, dfs)
}

as_logical_or_na <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  if (is.numeric(x) || is.integer(x)) {
    return(ifelse(is.na(x), NA, x != 0))
  }

  value <- trimws(tolower(as.character(x)))
  out <- rep(NA, length(value))
  out[value %in% c("true", "t", "1", "yes", "y", "nmd sensitive")] <- TRUE
  out[value %in% c("false", "f", "0", "no", "n", "nmd insensitive")] <- FALSE
  out
}

find_ptc_column <- function(df, required = TRUE) {
  candidates <- c(
    "PTC",
    "isoformFeatures__PTC",
    "orfAnalysis__PTC"
  )
  candidates <- unique(c(
    candidates[candidates %in% colnames(df)],
    grep("(^|__)PTC$", colnames(df), ignore.case = TRUE, value = TRUE)
  ))

  for (nm in candidates) {
    values <- as_logical_or_na(df[[nm]])
    if (any(!is.na(values))) {
      return(nm)
    }
  }

  if (isTRUE(required)) {
    stop(
      "No populated PTC column was found in the strict IsoformSwitchAnalyzeR object. ",
      "The cached object must contain ORF/PTC annotation before NMD status can be analyzed."
    )
  }
  NULL
}

resolve_strict_switch_cache <- function(fig6_analysis) {
  recorded <- fig6_analysis$meta$switch_rds %||% NA_character_
  if (
    length(recorded) == 1 &&
      !is.na(recorded) &&
      nzchar(recorded) &&
      file.exists(recorded)
  ) {
    resolved <- recorded
  } else {
    resolved <- find_isoform_switch_cache(prefer_strict = TRUE, required = TRUE)
  }

  if (identical(basename(resolved), "switchList_filtered_analyzed.rds")) {
    stop(
      "The Figure 6 analysis cache resolves to the non-strict IsoformSwitchAnalyzeR object: ",
      resolved,
      "\nRerun plot_fig_6.R with switchList_filtered_strict_analyzed.rds, then rerun ",
      "make_sup_tables.R so Figure 6 and the supplementary tables use the same calls."
    )
  }
  resolved
}

add_standard_comparison_columns <- function(df, reference_condition = "C1") {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(df) == 0) {
    return(df)
  }

  if ("condition_2" %in% colnames(df)) {
    df$timepoint <- as.character(df$condition_2)
  } else {
    comparison_col <- c("plotComparison", "Comparison")
    comparison_col <- comparison_col[comparison_col %in% colnames(df)]

    if (length(comparison_col) > 0) {
      comparison_text <- gsub(
        "[\r\n]+",
        " ",
        as.character(df[[comparison_col[[1]]]])
      )

      df$timepoint <- sub(
        paste0("^", reference_condition, "\\s+vs\\s+"),
        "",
        comparison_text
      )
    }
  }

  if (!"condition_1" %in% colnames(df)) {
    df$condition_1 <- reference_condition
  }

  if (!"condition_2" %in% colnames(df) && "timepoint" %in% colnames(df)) {
    df$condition_2 <- df$timepoint
  }

  if (!"comparison" %in% colnames(df) && "timepoint" %in% colnames(df)) {
    df$comparison <- paste0(df$timepoint, "_vs_", reference_condition)
  }

  df
}

extract_isoformswitch_nmd <- function(
  switch_list,
  annotation,
  alpha = 0.05,
  dIFcutoff = 0.10,
  reference_condition = "C1",
  timepoints = c("H1", "H3", "H24")
) {
  if (!requireNamespace("IsoformSwitchAnalyzeR", quietly = TRUE)) {
    stop("Package 'IsoformSwitchAnalyzeR' is required to extract NMD consequences.")
  }
  if (is.null(switch_list$isoformFeatures)) {
    stop("The strict IsoformSwitchAnalyzeR object does not contain $isoformFeatures.")
  }

  features <- .fig6_standardize_features(switch_list$isoformFeatures)
  keep <- features$condition_1 == reference_condition &
    features$condition_2 %in% timepoints
  if (!any(keep)) {
    stop("No ", reference_condition, "-to-hypoxia comparisons were found in the strict switch object.")
  }

  sub_obj <- .fig6_subset_switch_object(switch_list, keep)
  ptc_col <- find_ptc_column(sub_obj$isoformFeatures, required = TRUE)

  # This is a lightweight consequence calculation from the existing strict
  # switch calls and existing PTC annotation. It does not rerun switch detection.
  # Recomputing NMD alone with the same alpha and dIF cutoffs used by Figure 6
  # keeps the consequence universe aligned with the reported switch calls.
  nmd_obj <- IsoformSwitchAnalyzeR::analyzeSwitchConsequences(
    switchAnalyzeRlist = sub_obj,
    consequencesToAnalyze = "NMD_status",
    alpha = alpha,
    dIFcutoff = dIFcutoff,
    onlySigIsoforms = FALSE,
    removeNonConseqSwitches = FALSE,
    showProgress = FALSE,
    quiet = TRUE
  )

  details <- as.data.frame(
    nmd_obj$switchConsequence,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  details <- add_missing_columns(
    details,
    c(
      "gene_ref",
      "gene_id",
      "gene_name",
      "condition_1",
      "condition_2",
      "isoformUpregulated",
      "isoformDownregulated",
      "iso_ref_up",
      "iso_ref_down",
      "featureCompared",
      "isoformsDifferent",
      "switchConsequence"
    )
  )
  details <- details[
    as.character(details$featureCompared) == "NMD_status", ,
    drop = FALSE
  ]

  nmd_features <- as.data.frame(
    nmd_obj$isoformFeatures,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  ptc_col <- find_ptc_column(nmd_features, required = TRUE)
  ptc_values <- as_logical_or_na(nmd_features[[ptc_col]])

  if (nrow(details) > 0) {
    up_idx <- match(details$iso_ref_up, nmd_features$iso_ref)
    down_idx <- match(details$iso_ref_down, nmd_features$iso_ref)
    details$upregulated_isoform_PTC <- ptc_values[up_idx]
    details$downregulated_isoform_PTC <- ptc_values[down_idx]
    details$NMD_status_evaluable <- !is.na(details$upregulated_isoform_PTC) &
      !is.na(details$downregulated_isoform_PTC)
    details$NMD_change <- ifelse(
      !details$NMD_status_evaluable,
      "not_evaluable",
      ifelse(
        details$upregulated_isoform_PTC & !details$downregulated_isoform_PTC,
        "NMD_sensitive_gain",
        ifelse(
          !details$upregulated_isoform_PTC & details$downregulated_isoform_PTC,
          "NMD_sensitive_loss",
          "no_NMD_status_change"
        )
      )
    )
    details$timepoint <- as.character(details$condition_2)
    details$comparison <- paste0(details$condition_2, "_vs_", details$condition_1)
    details$PTC_source_column <- ptc_col

    up_annotation <- annotate_transcript_biotypes(
      data.frame(
        isoform_id = as.character(details$isoformUpregulated),
        stringsAsFactors = FALSE
      ),
      annotation = annotation,
      transcript_col = "isoform_id"
    )
    down_annotation <- annotate_transcript_biotypes(
      data.frame(
        isoform_id = as.character(details$isoformDownregulated),
        stringsAsFactors = FALSE
      ),
      annotation = annotation,
      transcript_col = "isoform_id"
    )
    details$upregulated_transcript_biotype <- up_annotation$transcript_biotype
    details$downregulated_transcript_biotype <- down_annotation$transcript_biotype
  } else {
    details$upregulated_isoform_PTC <- logical(0)
    details$downregulated_isoform_PTC <- logical(0)
    details$NMD_status_evaluable <- logical(0)
    details$NMD_change <- character(0)
    details$timepoint <- character(0)
    details$comparison <- character(0)
    details$PTC_source_column <- character(0)
    details$upregulated_transcript_biotype <- character(0)
    details$downregulated_transcript_biotype <- character(0)
  }

  consequence_summary <- IsoformSwitchAnalyzeR::extractConsequenceSummary(
    nmd_obj,
    consequencesToAnalyze = "NMD_status",
    includeCombined = FALSE,
    asFractionTotal = TRUE,
    alpha = alpha,
    dIFcutoff = dIFcutoff,
    plot = FALSE,
    plotGenes = TRUE,
    removeEmptyConsequences = FALSE,
    returnResult = TRUE
  )
  consequence_summary <- add_standard_comparison_columns(
    consequence_summary,
    reference_condition = reference_condition
  )
  consequence_summary$record_type <- rep("consequence_summary", nrow(consequence_summary))
  consequence_summary$count_unit <- rep("genes", nrow(consequence_summary))

  consequence_enrichment <- IsoformSwitchAnalyzeR::extractConsequenceEnrichment(
    nmd_obj,
    consequencesToAnalyze = "NMD_status",
    alpha = alpha,
    dIFcutoff = dIFcutoff,
    countGenes = TRUE,
    analysisOppositeConsequence = TRUE,
    plot = FALSE,
    returnResult = TRUE,
    returnSummary = TRUE
  )
  consequence_enrichment <- add_standard_comparison_columns(
    consequence_enrichment,
    reference_condition = reference_condition
  )
  consequence_enrichment <- add_missing_columns(
    consequence_enrichment,
    c("propQval", "propOfRelevantEvents")
  )
  consequence_enrichment$record_type <- rep(
    "consequence_enrichment",
    nrow(consequence_enrichment)
  )
  consequence_enrichment$count_unit <- rep("genes", nrow(consequence_enrichment))
  consequence_enrichment$selected_FDR <- rep(alpha, nrow(consequence_enrichment))
  consequence_enrichment$enriched_toward_NMD_sensitivity_at_selected_FDR <-
    !is.na(consequence_enrichment$propQval) &
      consequence_enrichment$propQval < alpha &
      consequence_enrichment$propOfRelevantEvents > 0.5
  consequence_enrichment$enriched_toward_NMD_insensitivity_at_selected_FDR <-
    !is.na(consequence_enrichment$propQval) &
      consequence_enrichment$propQval < alpha &
      consequence_enrichment$propOfRelevantEvents < 0.5

  pair_summary <- do.call(rbind, lapply(timepoints, function(tp) {
    x <- details[details$condition_2 == tp, , drop = FALSE]
    evaluable <- !is.na(x$NMD_status_evaluable) & x$NMD_status_evaluable
    changed <- x$NMD_change %in% c("NMD_sensitive_gain", "NMD_sensitive_loss")
    data.frame(
      record_type = "pairwise_PTC_coverage",
      count_unit = "pairwise_isoform_comparisons",
      timepoint = tp,
      comparison = paste0(tp, "_vs_", reference_condition),
      condition_1 = reference_condition,
      condition_2 = tp,
      n_pairwise_isoform_comparisons = nrow(x),
      n_PTC_evaluable_pairs = sum(evaluable),
      PTC_evaluable_fraction = if (nrow(x) > 0) sum(evaluable) / nrow(x) else NA_real_,
      n_NMD_status_changes = sum(changed, na.rm = TRUE),
      n_NMD_sensitive_gains = sum(x$NMD_change == "NMD_sensitive_gain", na.rm = TRUE),
      n_NMD_sensitive_losses = sum(x$NMD_change == "NMD_sensitive_loss", na.rm = TRUE),
      fraction_NMD_sensitive_gain_among_changes = if (sum(changed, na.rm = TRUE) > 0) {
        sum(x$NMD_change == "NMD_sensitive_gain", na.rm = TRUE) /
          sum(changed, na.rm = TRUE)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }))

  list(
    details = details,
    stats = rbind_fill_base(list(
      pair_summary,
      consequence_summary,
      consequence_enrichment
    )),
    analyzed_object = nmd_obj,
    ptc_source_column = ptc_col,
    alpha = alpha,
    dIFcutoff = dIFcutoff
  )
}

make_ptc_background_enrichment <- function(
  significant,
  background_features,
  timepoints = c("H1", "H3", "H24"),
  reference_condition = "C1",
  alpha = 0.10
) {
  significant <- as.data.frame(significant, stringsAsFactors = FALSE, check.names = FALSE)
  background_features <- .fig6_standardize_features(background_features)

  sig_ptc_col <- find_ptc_column(significant, required = TRUE)
  bg_ptc_col <- find_ptc_column(background_features, required = TRUE)
  significant$NMD_sensitive_by_PTC <- as_logical_or_na(significant[[sig_ptc_col]])
  background_features$NMD_sensitive_by_PTC <- as_logical_or_na(
    background_features[[bg_ptc_col]]
  )

  rows <- list()
  for (tp in timepoints) {
    sig_tp <- significant[
      significant$condition_1 == reference_condition &
        significant$condition_2 == tp, ,
      drop = FALSE
    ]
    bg_tp <- background_features[
      background_features$condition_1 == reference_condition &
        background_features$condition_2 == tp, ,
      drop = FALSE
    ]
    sig_tp <- sig_tp[!duplicated(sig_tp$isoform_id), , drop = FALSE]
    bg_tp <- bg_tp[!duplicated(bg_tp$isoform_id), , drop = FALSE]
    bg_eval <- bg_tp[!is.na(bg_tp$NMD_sensitive_by_PTC), , drop = FALSE]

    direction_sets <- list(
      all_switches = sig_tp,
      increased_in_hypoxia = sig_tp[
        sig_tp$usage_change_in_condition_2 == "increased" &
          !is.na(sig_tp$usage_change_in_condition_2), ,
        drop = FALSE
      ],
      decreased_in_hypoxia = sig_tp[
        sig_tp$usage_change_in_condition_2 == "decreased" &
          !is.na(sig_tp$usage_change_in_condition_2), ,
        drop = FALSE
      ]
    )

    rows[[tp]] <- do.call(rbind, lapply(names(direction_sets), function(direction_name) {
      sig_set <- direction_sets[[direction_name]]
      sig_eval <- sig_set[!is.na(sig_set$NMD_sensitive_by_PTC), , drop = FALSE]
      sig_eval <- sig_eval[
        sig_eval$isoform_id %in% bg_eval$isoform_id, ,
        drop = FALSE
      ]

      a <- sum(sig_eval$NMD_sensitive_by_PTC)
      b <- nrow(sig_eval) - a
      bg_positive <- sum(bg_eval$NMD_sensitive_by_PTC)
      c <- max(bg_positive - a, 0)
      d <- max((nrow(bg_eval) - bg_positive) - b, 0)
      ft <- if (nrow(sig_eval) > 0 && nrow(bg_eval) > nrow(sig_eval)) {
        stats::fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
      } else {
        NULL
      }

      data.frame(
        timepoint = tp,
        comparison = paste0(tp, "_vs_", reference_condition),
        switch_direction = direction_name,
        n_switching_transcripts = nrow(sig_set),
        n_switching_transcripts_with_PTC_status = nrow(sig_eval),
        switching_PTC_mapping_fraction = if (nrow(sig_set) > 0) {
          nrow(sig_eval) / nrow(sig_set)
        } else {
          NA_real_
        },
        n_NMD_sensitive_switching_transcripts = a,
        fraction_NMD_sensitive_switching = if (nrow(sig_eval) > 0) {
          a / nrow(sig_eval)
        } else {
          NA_real_
        },
        n_background_transcripts = nrow(bg_tp),
        n_background_transcripts_with_PTC_status = nrow(bg_eval),
        background_PTC_mapping_fraction = if (nrow(bg_tp) > 0) {
          nrow(bg_eval) / nrow(bg_tp)
        } else {
          NA_real_
        },
        n_NMD_sensitive_background_transcripts = bg_positive,
        fraction_NMD_sensitive_background = if (nrow(bg_eval) > 0) {
          bg_positive / nrow(bg_eval)
        } else {
          NA_real_
        },
        odds_ratio = if (is.null(ft)) NA_real_ else unname(ft$estimate),
        pvalue = if (is.null(ft)) NA_real_ else ft$p.value,
        significant_PTC_source_column = sig_ptc_col,
        background_PTC_source_column = bg_ptc_col,
        background_definition = paste(
          "All unique PTC-evaluable transcripts in the strict",
          "IsoformSwitchAnalyzeR isoformFeatures universe for this contrast"
        ),
        stringsAsFactors = FALSE
      )
    }))
  }

  out <- do.call(rbind, rows)
  out$padj_BH_all_timepoint_direction_tests <- stats::p.adjust(out$pvalue, method = "BH")
  out$enriched_for_NMD_sensitivity_at_fdr_0.10 <-
    !is.na(out$padj_BH_all_timepoint_direction_tests) &
      out$padj_BH_all_timepoint_direction_tests < alpha &
      out$odds_ratio > 1
  out$depleted_for_NMD_sensitivity_at_fdr_0.10 <-
    !is.na(out$padj_BH_all_timepoint_direction_tests) &
      out$padj_BH_all_timepoint_direction_tests < alpha &
      out$odds_ratio < 1
  rownames(out) <- NULL
  out
}

make_transcript_biotype_enrichment <- function(
  significant_by_tp,
  background_by_tp,
  annotation,
  method,
  timepoints = c("H1", "H3", "H24")
) {
  rows <- lapply(timepoints, function(tp) {
    sig <- as.data.frame(significant_by_tp[[tp]], stringsAsFactors = FALSE, check.names = FALSE)
    bg <- as.data.frame(background_by_tp[[tp]], stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(bg) == 0) {
      return(NULL)
    }
    if (ncol(sig) == 0) sig <- bg[0, , drop = FALSE]

    sig_tx_col <- .fig6_first_col(
      sig,
      c("feature_id", "enstx", "transcript_id", "isoform_id"),
      required = TRUE,
      label = paste(method, tp, "significant transcript ID")
    )
    bg_tx_col <- .fig6_first_col(
      bg,
      c("feature_id", "enstx", "transcript_id", "isoform_id"),
      required = TRUE,
      label = paste(method, tp, "background transcript ID")
    )

    sig$.transcript_key <- strip_ens_version(sig[[sig_tx_col]])
    bg$.transcript_key <- strip_ens_version(bg[[bg_tx_col]])
    sig <- sig[!duplicated(sig$.transcript_key), , drop = FALSE]
    bg <- bg[!duplicated(bg$.transcript_key), , drop = FALSE]
    sig <- annotate_transcript_biotypes(sig, annotation, transcript_col = sig_tx_col)
    bg <- annotate_transcript_biotypes(bg, annotation, transcript_col = bg_tx_col)

    n_sig_total <- nrow(sig)
    n_bg_total <- nrow(bg)
    sig_mapped <- sig[!is.na(sig$transcript_biotype), , drop = FALSE]
    bg_mapped <- bg[!is.na(bg$transcript_biotype), , drop = FALSE]
    if (nrow(bg_mapped) == 0) {
      return(NULL)
    }

    biotypes <- sort(unique(bg_mapped$transcript_biotype))
    out <- lapply(biotypes, function(bt) {
      sig_in <- sum(sig_mapped$transcript_biotype == bt)
      sig_n <- nrow(sig_mapped)
      bg_in <- sum(bg_mapped$transcript_biotype == bt)
      bg_n <- nrow(bg_mapped)
      a <- sig_in
      b <- sig_n - sig_in
      c <- max(bg_in - sig_in, 0)
      d <- max((bg_n - bg_in) - b, 0)
      ft <- if (sig_n > 0) {
        stats::fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
      } else {
        NULL
      }

      data.frame(
        method = method,
        timepoint = tp,
        comparison = paste0(tp, "_vs_C1"),
        transcript_biotype = bt,
        n_significant_transcripts = sig_in,
        n_significant_transcripts_with_biotype = sig_n,
        n_total_significant_transcripts = n_sig_total,
        significant_biotype_mapping_fraction = if (n_sig_total > 0) sig_n / n_sig_total else NA_real_,
        fraction_significant = if (sig_n > 0) sig_in / sig_n else NA_real_,
        n_background_transcripts = bg_in,
        n_background_transcripts_with_biotype = bg_n,
        n_total_background_transcripts = n_bg_total,
        background_biotype_mapping_fraction = bg_n / n_bg_total,
        fraction_background = bg_in / bg_n,
        odds_ratio = if (is.null(ft)) NA_real_ else unname(ft$estimate),
        pvalue = if (is.null(ft)) NA_real_ else ft$p.value,
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, out)
    out$padj <- stats::p.adjust(out$pvalue, method = "BH")
    out$enriched_at_fdr_0.10 <- !is.na(out$padj) & out$padj < 0.10 & out$odds_ratio > 1
    out$depleted_at_fdr_0.10 <- !is.na(out$padj) & out$padj < 0.10 & out$odds_ratio < 1
    out
  })
  rbind_fill_base(rows)
}


get_gene_id_col <- function(df) {
  candidates <- c("ensgene", "gene_id", "gene_id_full", "groupID")
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    return(NULL)
  }
  hit[[1]]
}

get_gene_vector <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(character(0))
  }
  col <- get_gene_id_col(df)
  if (is.null(col)) {
    return(character(0))
  }
  unique(strip_ens_version(df[[col]]))
}

get_symbol_map <- function(...) {
  dfs <- list(...)
  maps <- lapply(dfs, function(df) {
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }
    gene_col <- get_gene_id_col(df)
    if (is.null(gene_col) || !"symbol" %in% colnames(df)) {
      return(NULL)
    }
    out <- data.frame(
      ensgene = strip_ens_version(df[[gene_col]]),
      symbol = as.character(df$symbol),
      stringsAsFactors = FALSE
    )
    out <- out[!is.na(out$ensgene) & out$ensgene != "", , drop = FALSE]
    out <- out[!duplicated(out$ensgene), , drop = FALSE]
    out
  })
  maps <- Filter(Negate(is.null), maps)

  if (exists(".get_annot", mode = "function")) {
    annot <- tryCatch(.get_annot(), error = function(e) NULL)
    if (!is.null(annot) && all(c("ensgene", "symbol") %in% colnames(annot))) {
      maps <- c(list(data.frame(
        ensgene = strip_ens_version(annot$ensgene),
        symbol = as.character(annot$symbol),
        stringsAsFactors = FALSE
      )), maps)
    }
  }

  if (length(maps) == 0) {
    return(data.frame(ensgene = character(), symbol = character()))
  }
  out <- do.call(rbind, maps)
  out <- out[!is.na(out$ensgene) & out$ensgene != "", , drop = FALSE]
  out <- out[!duplicated(out$ensgene), , drop = FALSE]
  out
}

attach_symbol <- function(df, symbol_map) {
  if (is.null(df) || nrow(df) == 0) {
    return(df)
  }
  gene_col <- get_gene_id_col(df)
  if (is.null(gene_col)) {
    return(df)
  }
  if (!"ensgene" %in% colnames(df)) df$ensgene <- strip_ens_version(df[[gene_col]])
  if (!"symbol" %in% colnames(df) || all(is.na(df$symbol) | df$symbol == "")) {
    df$symbol <- symbol_map$symbol[match(strip_ens_version(df$ensgene), symbol_map$ensgene)]
  }
  df
}

make_sig_by_tp <- function(result_obj, slot = "gene_sig_by_tp") {
  x <- result_obj$results[[slot]]
  if (is.null(x)) {
    return(setNames(vector("list", length(timepoints)), timepoints))
  }
  out <- lapply(timepoints, function(tp) get_gene_vector(x[[tp]]))
  names(out) <- timepoints
  out
}


# ------------------------------------------------------------------------------
# Figure 1 GO enrichment helper
# ------------------------------------------------------------------------------
# The Figure 1 enrichment objects already contain one consolidated, labeled table
# at `enrichment$complete_tbl`. Export that table directly rather than recursively
# traversing `by_cluster`, `result`, and `raw$result`, which contain redundant
# copies of the same PANTHER results.

extract_fig1_go_complete_tbl <- function(fig1_go = NULL, fig1 = NULL) {
  candidates <- list(
    # Standalone Figure 1 enrichment object.
    fig1_go_enrichment_complete = if (!is.null(fig1_go)) {
      fig1_go$enrichment$complete_tbl
    } else {
      NULL
    },

    # Fallback if the standalone RDS stores complete_tbl at its top level.
    fig1_go_top_level_complete = if (!is.null(fig1_go)) {
      fig1_go$complete_tbl
    } else {
      NULL
    },

    # Figure 1 gene-level object produced by the updated plotting script.
    fig1_pathway_enrichment_complete = if (!is.null(fig1)) {
      fig1$pathway_enrichment$enrichment$complete_tbl
    } else {
      NULL
    },

    # Fallback for an alternate nested structure.
    fig1_enrichment_complete = if (!is.null(fig1)) {
      fig1$enrichment$complete_tbl
    } else {
      NULL
    }
  )

  for (nm in names(candidates)) {
    tbl <- candidates[[nm]]
    if (is.data.frame(tbl) && nrow(tbl) > 0) {
      message_step("Using Figure 1 GO table from: ", nm)
      tbl <- clean_for_excel(tbl)

      # Retain significant enrichment results at the cutoff used for Figure 1.
      adj_candidates <- names(tbl)[grepl("fdr|padj|adj.*p", names(tbl), ignore.case = TRUE)]
      if (length(adj_candidates) > 0) {
        adj <- suppressWarnings(as.numeric(tbl[[adj_candidates[[1]]]]))
        tbl <- tbl[is.na(adj) | adj < 0.10, , drop = FALSE]
      }

      # Keep a stable, reader-friendly order without changing the saved values.
      if ("cluster" %in% names(tbl)) {
        cluster_num <- suppressWarnings(as.numeric(sub(".*?([0-9]+)$", "\\1", as.character(tbl$cluster))))
        ontology_ord <- if ("ontology" %in% names(tbl)) match(tbl$ontology, c("BP", "MF")) else rep(NA_integer_, nrow(tbl))
        fdr_ord <- if (length(adj_candidates) > 0) suppressWarnings(as.numeric(tbl[[adj_candidates[[1]]]])) else rep(NA_real_, nrow(tbl))
        tbl <- tbl[order(cluster_num, ontology_ord, fdr_ord, na.last = TRUE), , drop = FALSE]
      }

      rownames(tbl) <- NULL
      return(tbl)
    }
  }

  NULL
}

# ------------------------------------------------------------------------------
# Locate result objects
# ------------------------------------------------------------------------------
analysis_dir <- file.path(PROJECT_ROOT, "results", "analysis")
plots_dir <- file.path(PROJECT_ROOT, "results", "plots")

paths <- list(
  deseq = c(
    file.path(analysis_dir, "deseq", "deseq_results.C1.rds"),
    Sys.glob(file.path(analysis_dir, "deseq", "deseq_results*.rds"))
  ),
  dexseq = c(
    file.path(analysis_dir, "dexseq", "dexseq_results.C1.summarizeOverlaps.multiOverlap.all.rds"),
    file.path(analysis_dir, "dexseq", "dexseq_results.C1.summarizeOverlaps.multiOverlap.unique.rds"),
    Sys.glob(file.path(analysis_dir, "dexseq", "dexseq_results.C1*.rds")),
    Sys.glob(file.path(analysis_dir, "dexseq", "dexseq_results*.rds"))
  ),
  drimseq = c(
    file.path(analysis_dir, "drimseq", "drimseq_results.C1.rds"),
    Sys.glob(file.path(analysis_dir, "drimseq", "drimseq_results*.rds"))
  ),
  suppa = c(
    file.path(analysis_dir, "suppa", "standard", "C1", "suppa_results.rds"),
    Sys.glob(file.path(analysis_dir, "suppa", "standard", "C1", "suppa_results*.rds")),
    Sys.glob(file.path(analysis_dir, "suppa", "**", "suppa_results*.rds"))
  ),
  fig1 = c(
    file.path(plots_dir, "fig_1", "fig_1_gene_level_analysis.rds"),
    Sys.glob(file.path(plots_dir, "fig_1", "fig_1_gene_level_analysis*.rds"))
  ),
  fig1_go = c(
    file.path(plots_dir, "fig_1", "fig_1_cluster_pathway_enrichment_analysis.rds"),
    file.path(plots_dir, "fig_1", "fig_1_cluster_GO_enrichment.rds"),
    Sys.glob(file.path(plots_dir, "fig_1", "*cluster*enrich*.rds")),
    Sys.glob(file.path(plots_dir, "fig_1", "*pathway*.rds")),
    Sys.glob(file.path(plots_dir, "fig_1", "*GO*.rds")),
    Sys.glob(file.path(plots_dir, "fig_1", "*go*.rds"))
  ),
  fig4_overlap = c(
    file.path(plots_dir, "fig_4", "fig_4_gene_vs_transcript_overlap_analysis.rds"),
    Sys.glob(file.path(plots_dir, "fig_4", "*overlap*.rds"))
  ),
  fig5_representation = c(
    file.path(plots_dir, "fig_5", "fig_5_splicing_pathway_representation_analysis.rds"),
    Sys.glob(file.path(plots_dir, "fig_5", "*representation*.rds"))
  ),
  fig5_enrichment = c(
    file.path(plots_dir, "fig_5", "fig_5_splicing_pathway_enrichment_analysis.rds"),
    Sys.glob(file.path(plots_dir, "fig_5", "*enrichment*.rds"))
  ),
  fig6_analysis = c(
    file.path(plots_dir, "fig_6", "fig_6_isoform_switch_analysis.rds"),
    Sys.glob(file.path(plots_dir, "fig_6", "*isoform*switch*analysis*.rds"))
  )
)

# Remove duplicate paths while preserving order.
paths <- lapply(paths, function(x) unique(x[nzchar(x)]))

deseq <- load_rds_optional(paths$deseq, "DESeq2 results")
dexseq <- load_rds_optional(paths$dexseq, "DEXSeq results")
drimseq <- load_rds_optional(paths$drimseq, "DRIMSeq results")
suppa <- load_rds_optional(paths$suppa, "SUPPA2 results")
fig1 <- load_rds_optional(paths$fig1, "Figure 1 gene-level cluster analysis")
fig1_go <- load_rds_optional(paths$fig1_go, "Figure 1 complete cluster-level GO enrichment")
fig4_overlap <- load_rds_optional(paths$fig4_overlap, "Figure 4 overlap analysis")
fig5_representation <- load_rds_optional(paths$fig5_representation, "Figure 5 pathway representation analysis")
fig5_enrichment <- load_rds_optional(paths$fig5_enrichment, "Figure 5 pathway enrichment analysis")
fig6_analysis <- load_rds_optional(paths$fig6_analysis, "Figure 6 isoform-switch analysis")

if (is.null(deseq) || is.null(dexseq) || is.null(drimseq) || is.null(suppa)) {
  stop("DESeq2, DEXSeq, DRIMSeq, and SUPPA2 result objects are required to build the core supplemental workbook.")
}

symbol_map <- get_symbol_map(
  deseq$results$gene_sig_all,
  dexseq$results$gene_sig_all,
  drimseq$results$gene_sig_all,
  suppa$results$gene_sig_all
)

# Use the exact project GTF, rather than a separately downloaded annotation, so
# transcript IDs and biotypes match the Salmon/DRIMSeq/IsoformSwitchAnalyzeR
# annotation universe.
message_step("Loading project GENCODE transcript and gene biotypes")
gencode_biotypes <- get_gencode_biotype_annotation()

# ------------------------------------------------------------------------------
# Build curated method-specific tables
# ------------------------------------------------------------------------------
message_step("Curating method-specific significant result tables")

deseq_degs <- attach_symbol(deseq$results$gene_sig_all, symbol_map)
deseq_degs <- add_missing_columns(deseq_degs, c("ensgene", "gene_id", "symbol", "timepoint", "comparison", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "direction"))
deseq_degs <- select_existing(deseq_degs, c("timepoint", "comparison", "ensgene", "gene_id", "symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "direction", "description", "biotype", "tool"))

dexseq_genes <- attach_symbol(dexseq$results$gene_sig_all, symbol_map)
dexseq_genes <- add_missing_columns(dexseq_genes, c("timepoint", "comparison", "ensgene", "gene_id", "symbol", "padj", "neglog10_padj", "significant"))
dexseq_genes <- select_existing(dexseq_genes, c("timepoint", "comparison", "ensgene", "gene_id", "gene_id_full", "symbol", "padj", "neglog10_padj", "significant", "description", "biotype", "tool"))

dexseq_exons <- attach_symbol(dexseq$results$exon_sig_all, symbol_map)
dexseq_exons <- add_missing_columns(dexseq_exons, c("timepoint", "comparison", "ensgene", "gene_id", "symbol", "exon_id", "feature_id", "log2fold", "pvalue", "padj", "significant"))
dexseq_exons <- select_existing(dexseq_exons, c("timepoint", "comparison", "ensgene", "gene_id", "gene_id_full", "symbol", "exon_id", "feature_id", "log2fold", "pvalue", "padj", "neglog10_padj", "significant", "description", "biotype", "tool"))

drimseq_genes <- attach_symbol(drimseq$results$gene_sig_all, symbol_map)
drimseq_genes <- add_missing_columns(drimseq_genes, c("timepoint", "comparison", "ensgene", "gene_id", "symbol", "pvalue", "adj_pvalue", "padj", "stageR_sig", "significant"))
drimseq_genes <- annotate_gene_biotypes(drimseq_genes, gencode_biotypes)
drimseq_genes <- select_existing(drimseq_genes, c("timepoint", "comparison", "ensgene", "gene_id", "gene_id_full", "gencode_gene_id", "symbol", "gene_symbol", "gene_biotype", "pvalue", "adj_pvalue", "padj", "neglog10_padj", "stageR_sig", "significant", "description", "biotype", "tool"))

drimseq_tx <- attach_symbol(drimseq$results$tx_sig_all, symbol_map)
drimseq_tx <- add_missing_columns(drimseq_tx, c("timepoint", "comparison", "ensgene", "gene_id", "symbol", "feature_id", "enstx", "pvalue", "adj_pvalue", "padj", "stageR_sig", "significant"))
drimseq_tx <- annotate_transcript_biotypes(
  drimseq_tx,
  annotation = gencode_biotypes,
  transcript_col = if ("feature_id" %in% colnames(drimseq_tx)) "feature_id" else "enstx"
)
drimseq_tx <- select_existing(drimseq_tx, c("timepoint", "comparison", "ensgene", "gene_id", "gene_id_full", "gencode_gene_id", "symbol", "gene_symbol", "feature_id", "enstx", "gencode_transcript_id", "transcript_id_full", "transcript_name", "transcript_biotype", "gene_biotype", "transcript_support_level", "pvalue", "adj_pvalue", "padj", "neglog10_padj", "stageR_sig", "significant", "description", "biotype", "tool"))

drimseq_biotype_enrichment <- make_transcript_biotype_enrichment(
  significant_by_tp = drimseq$results$tx_sig_by_tp,
  background_by_tp = drimseq$results$tx_full_by_tp,
  annotation = gencode_biotypes,
  method = "DRIMSeq",
  timepoints = timepoints
)

suppa_genes <- attach_symbol(suppa$results$gene_sig_all, symbol_map)
suppa_genes <- add_missing_columns(suppa_genes, c("timepoint", "comparison", "ensgene", "gene_id", "symbol"))
suppa_genes <- select_existing(suppa_genes, c("timepoint", "comparison", "ensgene", "gene_id", "ensgene_raw", "symbol", "description", "biotype", "tool"))

suppa_events <- attach_symbol(suppa$results$event_sig_all, symbol_map)
suppa_events <- add_missing_columns(suppa_events, c("timepoint", "comparison", "event_id", "ensgene", "gene_id", "symbol", "event_type", "chr", "strand", "dpsi", "pvalue", "direction"))
suppa_events <- select_existing(suppa_events, c("timepoint", "comparison", "event_id", "ensgene", "gene_id", "ensgene_raw", "symbol", "event_type", "chr", "coord1", "coord2", "coord3", "coord4", "strand", "dpsi", "pvalue", "neglog10_pvalue", "direction", "description", "biotype", "tool"))

# ------------------------------------------------------------------------------
# Integrated gene status table
# ------------------------------------------------------------------------------
message_step("Building integrated gene status table")

deg_by_tp <- make_sig_by_tp(deseq, "gene_sig_by_tp")
dex_by_tp <- make_sig_by_tp(dexseq, "gene_sig_by_tp")
drim_by_tp <- make_sig_by_tp(drimseq, "gene_sig_by_tp")
suppa_by_tp <- make_sig_by_tp(suppa, "gene_sig_by_tp")
asg_by_tp <- lapply(timepoints, function(tp) unique(c(dex_by_tp[[tp]], drim_by_tp[[tp]], suppa_by_tp[[tp]])))
names(asg_by_tp) <- timepoints

all_genes <- unique(unlist(c(deg_by_tp, dex_by_tp, drim_by_tp, suppa_by_tp), use.names = FALSE))
all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]

integrated <- data.frame(ensgene = sort(all_genes), stringsAsFactors = FALSE)
integrated$symbol <- symbol_map$symbol[match(integrated$ensgene, symbol_map$ensgene)]

for (tp in timepoints) {
  integrated[[paste0(tp, "_DEG")]] <- integrated$ensgene %in% deg_by_tp[[tp]]
  integrated[[paste0(tp, "_DEXSeq")]] <- integrated$ensgene %in% dex_by_tp[[tp]]
  integrated[[paste0(tp, "_DRIMSeq")]] <- integrated$ensgene %in% drim_by_tp[[tp]]
  integrated[[paste0(tp, "_SUPPA2")]] <- integrated$ensgene %in% suppa_by_tp[[tp]]
  integrated[[paste0(tp, "_ASG_combined")]] <- integrated$ensgene %in% asg_by_tp[[tp]]
  integrated[[paste0(tp, "_classification")]] <- ifelse(
    integrated[[paste0(tp, "_DEG")]] & integrated[[paste0(tp, "_ASG_combined")]], "DEG_and_ASG",
    ifelse(integrated[[paste0(tp, "_DEG")]], "DEG_only",
      ifelse(integrated[[paste0(tp, "_ASG_combined")]], "ASG_only", "not_detected")
    )
  )
}

integrated$n_DEG_timepoints <- rowSums(integrated[paste0(timepoints, "_DEG")])
integrated$n_ASG_timepoints <- rowSums(integrated[paste0(timepoints, "_ASG_combined")])
integrated$any_DEG <- integrated$n_DEG_timepoints > 0
integrated$any_ASG <- integrated$n_ASG_timepoints > 0
integrated$detected_layers <- ifelse(integrated$any_DEG & integrated$any_ASG, "both_layers",
  ifelse(integrated$any_DEG, "gene_level_only", "transcript_level_only")
)

# ------------------------------------------------------------------------------
# Figure 4C transition table
# ------------------------------------------------------------------------------
message_step("Building regulatory transition table")

transitions <- data.frame(note = "Figure 4 overlap analysis object was not found; transition table not generated.")
if (!is.null(fig4_overlap) && !is.null(fig4_overlap$results$highlighted_sets)) {
  highlighted_sets <- fig4_overlap$results$highlighted_sets
  alluvial_df <- fig4_overlap$results$alluvial_df %||% data.frame()
  trans_list <- lapply(names(highlighted_sets), function(group_name) {
    genes <- strip_ens_version(highlighted_sets[[group_name]])
    if (length(genes) == 0) {
      return(NULL)
    }
    meta <- alluvial_df[alluvial_df$Group == group_name, , drop = FALSE]
    data.frame(
      group = group_name,
      ensgene = genes,
      symbol = symbol_map$symbol[match(genes, symbol_map$ensgene)],
      H1_state = if (nrow(meta) > 0) meta$H1[[1]] else NA_character_,
      H3_state = if (nrow(meta) > 0) meta$H3[[1]] else NA_character_,
      H24_state = if (nrow(meta) > 0) meta$H24[[1]] else NA_character_,
      group_size = if (nrow(meta) > 0) meta$Freq[[1]] else length(genes),
      stringsAsFactors = FALSE
    )
  })
  transitions <- do.call(rbind, Filter(Negate(is.null), trans_list))
  if (is.null(transitions)) transitions <- data.frame(note = "No highlighted transition genes found.")
}

# ------------------------------------------------------------------------------
# Figure 1 cluster table
# ------------------------------------------------------------------------------
message_step("Building Figure 1 cluster table")

fig1_clusters <- data.frame(note = "Figure 1 cluster analysis object was not found; cluster table not generated.")
if (!is.null(fig1) && !is.null(fig1$clustering$clusters_df)) {
  fig1_clusters <- fig1$clustering$clusters_df
  fig1_clusters$ensgene <- strip_ens_version(fig1_clusters$ensgene)
  fig1_clusters$symbol <- symbol_map$symbol[match(fig1_clusters$ensgene, symbol_map$ensgene)]

  expr <- fig1$expression$counts_normalized_average_scaled_sig
  if (!is.null(expr)) {
    expr_df <- as.data.frame(expr, stringsAsFactors = FALSE, check.names = FALSE)
    expr_df$ensgene <- strip_ens_version(rownames(expr_df))
    rownames(expr_df) <- NULL
    fig1_clusters <- merge(fig1_clusters, expr_df, by = "ensgene", all.x = TRUE)
  }
  fig1_clusters <- fig1_clusters[order(fig1_clusters$cluster, fig1_clusters$ensgene), , drop = FALSE]
}


# ------------------------------------------------------------------------------
# Complete Figure 1 cluster-level GO enrichment table
# ------------------------------------------------------------------------------
message_step("Recovering complete Figure 1 cluster-level GO enrichment results")

fig1_go_full <- extract_fig1_go_complete_tbl(fig1_go = fig1_go, fig1 = fig1)
if (is.null(fig1_go_full) || nrow(fig1_go_full) == 0) {
  stop(
    "Could not recover the complete Figure 1 GO enrichment table. ",
    "Save the full cluster-level enrichment result as an RDS file under ",
    file.path(plots_dir, "fig_1"),
    " (for example, fig_1_cluster_pathway_enrichment_analysis.rds). ",
    "The script intentionally stops rather than silently exporting only the 15 curated Figure 1C terms."
  )
}

if (nrow(fig1_go_full) <= 15) {
  stop(
    "Only ", nrow(fig1_go_full),
    " Figure 1 GO rows were recovered. This appears to be the curated Figure 1C subset, ",
    "not the complete cluster enrichment analysis. Please save/load the full enrichment object."
  )
}
message_step("Recovered ", nrow(fig1_go_full), " complete Figure 1 GO enrichment rows")

# ------------------------------------------------------------------------------
# Pathway tables
# ------------------------------------------------------------------------------
message_step("Building pathway tables")

pathway_representation <- data.frame(note = "Figure 5 pathway representation object was not found.")
if (!is.null(fig5_representation) && !is.null(fig5_representation$pathway_representation$wide)) {
  pathway_representation <- fig5_representation$pathway_representation$wide
}

pathway_asg <- data.frame(note = "Figure 5 ASG pathway enrichment object was not found.")
pathway_deg <- data.frame(note = "Figure 5 DEG pathway enrichment object was not found.")
pathway_overlap <- data.frame(note = "Figure 5 pathway overlap object was not found.")

if (!is.null(fig5_enrichment)) {
  if (!is.null(fig5_enrichment$enrichment$as_pathways_tbl) && nrow(fig5_enrichment$enrichment$as_pathways_tbl) > 0) {
    pathway_asg <- fig5_enrichment$enrichment$as_pathways_tbl
  }
  if (!is.null(fig5_enrichment$enrichment$deg_pathways_tbl) && nrow(fig5_enrichment$enrichment$deg_pathways_tbl) > 0) {
    pathway_deg <- fig5_enrichment$enrichment$deg_pathways_tbl
  }

  overlap_tables <- list()
  if (!is.null(fig5_enrichment$overlap$as_pathways_overlap_df)) {
    x <- fig5_enrichment$overlap$as_pathways_overlap_df
    x$table <- "ASG pathway overlap across timepoints"
    overlap_tables[[length(overlap_tables) + 1]] <- x
  }
  if (!is.null(fig5_enrichment$overlap$pathway_compare_df)) {
    x <- fig5_enrichment$overlap$pathway_compare_df
    x$table <- "ASG vs DEG pathway overlap by timepoint"
    overlap_tables[[length(overlap_tables) + 1]] <- x
  }
  if (length(overlap_tables) > 0) pathway_overlap <- rbind_fill_base(overlap_tables)
}

# ------------------------------------------------------------------------------
# Figure 6 IsoformSwitchAnalyzeR and transcript-biotype tables
# ------------------------------------------------------------------------------
message_step("Building IsoformSwitchAnalyzeR and transcript-biotype tables")

fig6_recorded_switch_rds <- if (is.null(fig6_analysis)) {
  NA_character_
} else {
  fig6_analysis$meta$switch_rds %||% NA_character_
}
fig6_known_nonstrict <- length(fig6_recorded_switch_rds) == 1 &&
  !is.na(fig6_recorded_switch_rds) &&
  identical(basename(fig6_recorded_switch_rds), "switchList_filtered_analyzed.rds")

if (
  is.null(fig6_analysis) ||
    !identical(fig6_analysis$meta$analysis_version, .fig6_analysis_version) ||
    is.null(fig6_analysis$tables$isoform_switches) ||
    is.null(fig6_analysis$tables$biotype_enrichment) ||
    fig6_known_nonstrict
) {
  switch_rds <- find_isoform_switch_cache(prefer_strict = TRUE, required = TRUE)
  if (identical(basename(switch_rds), "switchList_filtered_analyzed.rds")) {
    stop(
      "Could not locate switchList_filtered_strict_analyzed.rds. ",
      "The supplementary IsoformSwitchAnalyzeR results must match the strict Figure 6 analysis."
    )
  }
  fig6_analysis <- run_isoform_switch_analysis(
    switch_rds = switch_rds,
    deseq_results = deseq,
    outdir = file.path(plots_dir, "fig_6"),
    force = FALSE,
    save_tables = TRUE
  )
}

strict_switch_rds <- resolve_strict_switch_cache(fig6_analysis)
message_step("Loading strict IsoformSwitchAnalyzeR object for NMD/PTC analysis: ", strict_switch_rds)
strict_switch_list <- readRDS(strict_switch_rds)

isoform_switches <- as.data.frame(
  fig6_analysis$tables$isoform_switches,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (nrow(isoform_switches) == 0) {
  isoform_switches <- data.frame(note = "No significant isoform switches were recovered.")
} else {
  preferred_switch_cols <- c(
    "timepoint",
    "comparison",
    "condition_1",
    "condition_2",
    "gene_id",
    "gene_name",
    "gene_symbol",
    "isoform_id",
    "gencode_transcript_id",
    "transcript_id_full",
    "transcript_name",
    "transcript_biotype",
    "gene_biotype",
    "transcript_support_level",
    "dIF",
    "usage_change_in_condition_2",
    "gene_log2_fold_change",
    "DEG",
    "iso_q_value",
    "gene_q_value",
    "isoform_switch_q_value",
    "isoform_switch_P_value",
    "gene_switch_q_value",
    "gencode_NMD_biotype"
  )
  nmd_or_orf_cols <- grep(
    "nmd|ptc|premature|orf|coding|cpc|consequence",
    colnames(isoform_switches),
    ignore.case = TRUE,
    value = TRUE
  )
  isoform_ptc_col <- find_ptc_column(isoform_switches, required = TRUE)
  isoform_switches$NMD_sensitive_by_PTC <- as_logical_or_na(
    isoform_switches[[isoform_ptc_col]]
  )
  isoform_switches$PTC_source_column <- isoform_ptc_col
  isoform_switches <- select_existing(
    isoform_switches,
    unique(c(
      preferred_switch_cols,
      "NMD_sensitive_by_PTC",
      "PTC_source_column",
      nmd_or_orf_cols
    ))
  )
}

isoform_biotype_enrichment <- as.data.frame(
  fig6_analysis$tables$biotype_enrichment,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (nrow(isoform_biotype_enrichment) == 0) {
  isoform_biotype_enrichment <- data.frame(
    note = paste(
      "No transcript-biotype enrichment results were recovered.",
      "The background should be all eligible transcripts in isoformFeatures for each contrast."
    )
  )
}

isoform_switch_summary <- fig6_analysis$summary
if (is.null(isoform_switch_summary) || !is.data.frame(isoform_switch_summary)) {
  isoform_switch_summary <- data.frame()
}

message_step("Extracting strict-cache NMD switch consequences and enrichment")
nmd_results <- extract_isoformswitch_nmd(
  switch_list = strict_switch_list,
  annotation = gencode_biotypes,
  alpha = 0.05,
  dIFcutoff = 0.10,
  reference_condition = "C1",
  timepoints = timepoints
)
isoform_nmd_details <- nmd_results$details
isoform_nmd_stats <- nmd_results$stats

message_step("Testing PTC/NMD-sensitivity enrichment against the expressed transcript background")
isoform_ptc_summary <- make_ptc_background_enrichment(
  significant = fig6_analysis$tables$isoform_switches,
  background_features = strict_switch_list$isoformFeatures,
  timepoints = timepoints,
  reference_condition = "C1",
  alpha = 0.10
)

# ------------------------------------------------------------------------------
# Workbook assembly
# ------------------------------------------------------------------------------
message_step("Writing Excel workbook")

dir.create(dirname(OUT_XLSX), recursive = TRUE, showWarnings = FALSE)
wb <- openxlsx::createWorkbook()

readme <- data.frame(
  field = c(
    "Manuscript title",
    "Workbook purpose",
    "Included results",
    "Biotype enrichment background",
    "NMD consequence analysis",
    "PTC background enrichment",
    "Strict IsoformSwitchAnalyzeR source",
    "Data availability note",
    "Thresholds",
    "Timepoints",
    "Generated on",
    "Project root"
  ),
  value = c(
    "Gene-level RNA-seq obscures extensive transcript-level remodeling during hypoxic adaptation",
    "Curated supplementary data tables supporting DEG, alternative splicing/transcript usage, isoform switching, GENCODE transcript-biotype enrichment, integrated gene classification, regulatory transition, and pathway analyses.",
    "Significant DESeq2 DEGs; significant DEXSeq genes/exons; significant DRIMSeq genes/transcripts annotated with GENCODE gene/transcript biotypes; significant SUPPA2 genes/events; strict-cache IsoformSwitchAnalyzeR calls, pairwise NMD consequences, NMD consequence enrichment, and PTC-background enrichment; DRIMSeq and IsoformSwitchAnalyzeR transcript-biotype enrichment; integrated gene status; Figure 4 transition groups; Figure 1 DEG clusters; complete significant Figure 1 cluster-level GO enrichment results; Figure 5 pathway summaries.",
    "For each C1-to-hypoxia contrast, the IsoformSwitchAnalyzeR background is all unique transcripts in isoformFeatures and the DRIMSeq background is all unique transcripts in tx_full_by_tp. Thus each test uses transcripts eligible for that method/contrast; GENCODE v45 supplies annotation labels but does not define the tested background.",
    "NMD status was recalculated in memory from the strict cached switch calls and their existing PTC annotation using analyzeSwitchConsequences(consequencesToAnalyze = 'NMD_status', alpha = 0.05, dIFcutoff = 0.10, onlySigIsoforms = FALSE), matching the Figure 6 extractTopSwitches defaults. NMD-sensitive gain means the isoform used more in hypoxia is PTC-positive while the paired isoform used less is PTC-negative. The package consequence-enrichment test counts genes and uses a binomial null proportion of 0.5, with BH correction reported by IsoformSwitchAnalyzeR.",
    "For each timepoint and switch direction, PTC-positive (predicted NMD-sensitive) switching transcripts were compared with all other PTC-evaluable transcripts in the strict isoformFeatures background using a two-sided Fisher exact test. BH correction is applied jointly across the timepoint-by-direction tests. This analysis is distinct from the GENCODE nonsense_mediated_decay biotype test.",
    strict_switch_rds,
    "Full unfiltered method-level statistical result universes remain in the analysis RDS objects. The workbook contains significant calls and the complete Figure 1 cluster-level GO results.",
    "Primary significance threshold used throughout manuscript: FDR or adjusted P < 0.10 unless otherwise noted by tool-specific output. IsoformSwitchAnalyzeR Figure 6 calls and NMD consequences use the package default FDR < 0.05 plus |dIF| > 0.10; the supplemental PTC-background enrichment uses BH-adjusted P < 0.10.",
    paste(paste(names(timepoint_labels), timepoint_labels, sep = " = "), collapse = "; "),
    as.character(Sys.time()),
    PROJECT_ROOT
  ),
  stringsAsFactors = FALSE
)

method_summary <- rbind_fill_base(list(
  cbind(method = "DESeq2", clean_for_excel(deseq$summary)),
  cbind(method = "DEXSeq", clean_for_excel(dexseq$summary)),
  cbind(method = "DRIMSeq", clean_for_excel(drimseq$summary)),
  cbind(method = "SUPPA2", clean_for_excel(suppa$summary)),
  cbind(method = "IsoformSwitchAnalyzeR", clean_for_excel(isoform_switch_summary))
))

add_sheet(wb, "01_README", readme)
add_sheet(wb, "02_Method_summary", method_summary)
add_sheet(wb, "03_DESeq2_DEGs", deseq_degs)
add_sheet(wb, "04_DEXSeq_genes", dexseq_genes)
add_sheet(wb, "05_DEXSeq_exons", dexseq_exons)
add_sheet(wb, "06_DRIMSeq_genes", drimseq_genes)
add_sheet(wb, "07_DRIMSeq_transcripts", drimseq_tx)
add_sheet(wb, "08_SUPPA2_genes", suppa_genes)
add_sheet(wb, "09_SUPPA2_events", suppa_events)
add_sheet(wb, "10_Integrated_gene_status", integrated)
add_sheet(wb, "11_Regulatory_transitions", transitions)
add_sheet(wb, "12_DEG_clusters_Fig1", fig1_clusters)
add_sheet(wb, "13_Fig1_GO_complete", fig1_go_full)
add_sheet(wb, "14_Pathway_representation", pathway_representation)
add_sheet(wb, "15_Pathway_enrichment_ASG", pathway_asg)
add_sheet(wb, "16_Pathway_enrichment_DEG", pathway_deg)
add_sheet(wb, "17_Pathway_overlap_summary", pathway_overlap)
add_sheet(wb, "18_IsoformSwitch_calls", isoform_switches)
add_sheet(wb, "19_IsoformSwitch_biotype", isoform_biotype_enrichment)
add_sheet(wb, "20_DRIMSeq_biotype", drimseq_biotype_enrichment)
add_sheet(wb, "21_IsoformSwitch_NMD_pairs", isoform_nmd_details)
add_sheet(wb, "22_IsoformSwitch_NMD_stats", isoform_nmd_stats)
add_sheet(wb, "23_IsoformSwitch_PTC_enrich", isoform_ptc_summary)

openxlsx::saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
message_step("Done: ", OUT_XLSX)
