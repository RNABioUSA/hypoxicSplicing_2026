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
#   $PROJECT_ROOT/results/salmon/<condition>/<sample>/quant.sf
#   $PROJECT_ROOT/results/star/**/Log.final.out
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
# Optional QC overrides:
#   export SALMON_COLDATA_TSV=/path/to/results/salmon/coldata.tsv
#   export SALMON_QC_DIR=/path/to/results/salmon
#   export STAR_QC_DIR=/path/to/results/star
#   export SAMPLE_METADATA_TSV=/path/to/config/metadata/samples.tsv
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
  if (exists(".strip_ens_version", mode = "function")) return(.strip_ens_version(x))
  sub("\\.[0-9]+$", "", as.character(x))
}

# ------------------------------------------------------------------------------
# General helpers
# ------------------------------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x

message_step <- function(...) message("[supplement] ", ...)

first_existing <- function(paths, required = FALSE, label = "file") {
  hit <- paths[file.exists(paths)]
  if (length(hit) > 0) return(hit[[1]])
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
  if (is.null(df)) return(data.frame())
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(df) == 0 && ncol(df) == 0) return(data.frame(note = "No records available."))

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
  integer_style <- openxlsx::createStyle(numFmt = "#,##0")
  percent_style <- openxlsx::createStyle(numFmt = "0.00")

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
    if (nrow(df) > 0) {
      numeric_cols <- vapply(df, is.numeric, logical(1))
      count_cols <- which(
        numeric_cols &
          grepl(
            "count|fragment|mapped.*number|targets|n_features|n_genes",
            names(df),
            ignore.case = TRUE
          )
      )
      percent_cols <- which(
        numeric_cols &
          grepl("(_pct|percent|percentage)", names(df), ignore.case = TRUE)
      )
      if (length(count_cols) > 0) {
        openxlsx::addStyle(
          wb,
          sheet_name,
          integer_style,
          rows = 2:(nrow(df) + 1),
          cols = count_cols,
          gridExpand = TRUE,
          stack = TRUE
        )
      }
      if (length(percent_cols) > 0) {
        openxlsx::addStyle(
          wb,
          sheet_name,
          percent_style,
          rows = 2:(nrow(df) + 1),
          cols = percent_cols,
          gridExpand = TRUE,
          stack = TRUE
        )
      }
    }
  }

  openxlsx::freezePane(wb, sheet_name, firstActiveRow = freeze_row + 1, firstActiveCol = freeze_col + 1)
  invisible(sheet_name)
}

select_existing <- function(df, cols) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  cols <- cols[cols %in% colnames(df)]
  if (length(cols) == 0) return(df)
  df[, unique(cols), drop = FALSE]
}

add_missing_columns <- function(df, cols) {
  for (nm in cols) if (!nm %in% colnames(df)) df[[nm]] <- NA
  df
}

first_present_column <- function(df, candidates, require_data = FALSE) {
  candidates <- candidates[candidates %in% colnames(df)]
  if (length(candidates) == 0) return(NULL)
  if (!isTRUE(require_data)) return(candidates[[1]])
  for (nm in candidates) {
    if (any(!is.na(df[[nm]]))) return(nm)
  }
  NULL
}

tested_row_mask <- function(df, statistic_candidates = c(
                              "padj", "adj_pvalue", "pvalue",
                              "iso_q_value", "isoform_switch_q_value"
                            )) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(logical())
  statistic_col <- first_present_column(df, statistic_candidates, require_data = TRUE)
  if (is.null(statistic_col)) return(rep(TRUE, nrow(df)))
  !is.na(df[[statistic_col]])
}

count_tested_rows <- function(df, statistic_candidates = c(
                                "padj", "adj_pvalue", "pvalue",
                                "iso_q_value", "isoform_switch_q_value"
                              )) {
  sum(tested_row_mask(df, statistic_candidates))
}

count_unique_genes <- function(
  df,
  tested_only = TRUE,
  statistic_candidates = c(
    "padj", "adj_pvalue", "pvalue",
    "iso_q_value", "isoform_switch_q_value"
  )
) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(0L)
  gene_col <- first_present_column(
    df,
    c("ensgene", "gene_id", "gene_id_full", "gencode_gene_id")
  )
  if (is.null(gene_col)) return(NA_integer_)
  keep <- if (isTRUE(tested_only)) {
    tested_row_mask(df, statistic_candidates = statistic_candidates)
  } else {
    rep(TRUE, nrow(df))
  }
  genes <- strip_ens_version(df[[gene_col]][keep])
  as.integer(length(unique(genes[!is.na(genes) & nzchar(genes)])))
}

as_numeric_or_na <- function(x) {
  out <- suppressWarnings(as.numeric(gsub("%", "", as.character(x), fixed = TRUE)))
  ifelse(is.finite(out), out, NA_real_)
}

normalize_sample_key <- function(x) {
  tolower(gsub("[^A-Za-z0-9]+", "", as.character(x)))
}

resolve_relative_path <- function(path, root = PROJECT_ROOT) {
  if (is.na(path) || !nzchar(path)) return(NA_character_)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) path <- file.path(root, path)
  normalizePath(path, mustWork = FALSE)
}

# rbind helper that preserves all columns across data frames.
rbind_fill_base <- function(dfs) {
  dfs <- Filter(function(x) !is.null(x) && is.data.frame(x), dfs)
  if (length(dfs) == 0) return(data.frame())
  cols <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  dfs <- lapply(dfs, function(x) {
    missing <- setdiff(cols, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[, cols, drop = FALSE]
  })
  do.call(rbind, dfs)
}

as_logical_or_na <- function(x) {
  if (is.logical(x)) return(x)
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
    if (any(!is.na(values))) return(nm)
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
  if (nrow(df) == 0) return(df)

  if ("condition_2" %in% colnames(df)) {
    df$timepoint <- as.character(df$condition_2)
  } else if ("plotComparison" %in% colnames(df)) {
    comparison_text <- gsub("[\r\n]+", " ", as.character(df$plotComparison))
    df$timepoint <- sub(
      paste0("^", reference_condition, "\\s+vs\\s+"),
      "",
      comparison_text
    )
  }
  if (!"condition_1" %in% colnames(df)) df$condition_1 <- reference_condition
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
    as.character(details$featureCompared) == "NMD_status",
    ,
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
    analysisOppositeConsequence = FALSE,
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
        significant$condition_2 == tp,
      ,
      drop = FALSE
    ]
    bg_tp <- background_features[
      background_features$condition_1 == reference_condition &
        background_features$condition_2 == tp,
      ,
      drop = FALSE
    ]
    sig_tp <- sig_tp[!duplicated(sig_tp$isoform_id), , drop = FALSE]
    bg_tp <- bg_tp[!duplicated(bg_tp$isoform_id), , drop = FALSE]
    bg_eval <- bg_tp[!is.na(bg_tp$NMD_sensitive_by_PTC), , drop = FALSE]

    direction_sets <- list(
      all_switches = sig_tp,
      increased_in_hypoxia = sig_tp[
        sig_tp$usage_change_in_condition_2 == "increased" &
          !is.na(sig_tp$usage_change_in_condition_2),
        ,
        drop = FALSE
      ],
      decreased_in_hypoxia = sig_tp[
        sig_tp$usage_change_in_condition_2 == "decreased" &
          !is.na(sig_tp$usage_change_in_condition_2),
        ,
        drop = FALSE
      ]
    )

    rows[[tp]] <- do.call(rbind, lapply(names(direction_sets), function(direction_name) {
      sig_set <- direction_sets[[direction_name]]
      sig_eval <- sig_set[!is.na(sig_set$NMD_sensitive_by_PTC), , drop = FALSE]
      sig_eval <- sig_eval[
        sig_eval$isoform_id %in% bg_eval$isoform_id,
        ,
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
    if (nrow(bg) == 0) return(NULL)
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
    if (nrow(bg_mapped) == 0) return(NULL)

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
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

get_gene_vector <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(character(0))
  col <- get_gene_id_col(df)
  if (is.null(col)) return(character(0))
  unique(strip_ens_version(df[[col]]))
}

get_symbol_map <- function(...) {
  dfs <- list(...)
  maps <- lapply(dfs, function(df) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    gene_col <- get_gene_id_col(df)
    if (is.null(gene_col) || !"symbol" %in% colnames(df)) return(NULL)
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

  if (length(maps) == 0) return(data.frame(ensgene = character(), symbol = character()))
  out <- do.call(rbind, maps)
  out <- out[!is.na(out$ensgene) & out$ensgene != "", , drop = FALSE]
  out <- out[!duplicated(out$ensgene), , drop = FALSE]
  out
}

attach_symbol <- function(df, symbol_map) {
  if (is.null(df) || nrow(df) == 0) return(df)
  gene_col <- get_gene_id_col(df)
  if (is.null(gene_col)) return(df)
  if (!"ensgene" %in% colnames(df)) df$ensgene <- strip_ens_version(df[[gene_col]])
  if (!"symbol" %in% colnames(df) || all(is.na(df$symbol) | df$symbol == "")) {
    df$symbol <- symbol_map$symbol[match(strip_ens_version(df$ensgene), symbol_map$ensgene)]
  }
  df
}

make_sig_by_tp <- function(result_obj, slot = "gene_sig_by_tp") {
  x <- result_obj$results[[slot]]
  if (is.null(x)) return(setNames(vector("list", length(timepoints)), timepoints))
  out <- lapply(timepoints, function(tp) get_gene_vector(x[[tp]]))
  names(out) <- timepoints
  out
}

# ------------------------------------------------------------------------------
# Standardized method-summary helpers
# ------------------------------------------------------------------------------
make_standard_method_summary <- function(
  deseq,
  dexseq,
  drimseq,
  suppa,
  fig6_analysis,
  strict_switch_list,
  timepoints = c("H1", "H3", "H24"),
  reference_condition = "C1"
) {
  get_tp_table <- function(object, slot, tp) {
    x <- object$results[[slot]]
    if (is.null(x) || is.null(x[[tp]])) return(data.frame())
    as.data.frame(x[[tp]], stringsAsFactors = FALSE, check.names = FALSE)
  }

  summary_row <- function(
    method,
    tp,
    features,
    significance_metric,
    significance_cutoff,
    effect_size_cutoff = NA_character_,
    n_features_tested = NA_integer_,
    n_sig_features = NA_integer_,
    n_genes_tested = NA_integer_,
    n_sig_genes = NA_integer_,
    notes = NA_character_
  ) {
    data.frame(
      method = method,
      timepoint = tp,
      comparison = paste(tp, "vs", reference_condition),
      features = features,
      significance_metric = significance_metric,
      significance_cutoff = as.character(significance_cutoff),
      effect_size_cutoff = as.character(effect_size_cutoff),
      n_features_tested = as.integer(n_features_tested),
      n_sig_features = as.integer(n_sig_features),
      n_genes_tested = as.integer(n_genes_tested),
      n_sig_genes = as.integer(n_sig_genes),
      notes = notes,
      stringsAsFactors = FALSE
    )
  }

  rows <- list()
  add_row <- function(x) rows[[length(rows) + 1L]] <<- x

  for (tp in timepoints) {
    # DESeq2: genes are both the tested features and the reported significant genes.
    full <- get_tp_table(deseq, "gene_full_by_tp", tp)
    sig <- get_tp_table(deseq, "gene_sig_by_tp", tp)
    tested_genes <- deseq$results$tested_universe_by_tp[[tp]] %||% character()
    n_tested <- if (length(tested_genes) > 0) {
      length(unique(strip_ens_version(tested_genes)))
    } else {
      count_tested_rows(full, c("padj", "pvalue"))
    }
    add_row(summary_row(
      method = "DESeq2",
      tp = tp,
      features = "genes",
      significance_metric = "BH-adjusted P value",
      significance_cutoff = deseq$meta$padj_cutoff %||% 0.10,
      n_features_tested = n_tested,
      n_sig_features = nrow(sig),
      n_genes_tested = n_tested,
      n_sig_genes = nrow(sig),
      notes = "Independent-filtered genes with non-missing adjusted P values define the tested universe."
    ))

    # DEXSeq: exon bins are features; genes summarize one or more significant bins.
    exon_full <- get_tp_table(dexseq, "exon_full_by_tp", tp)
    exon_sig <- get_tp_table(dexseq, "exon_sig_by_tp", tp)
    gene_full <- get_tp_table(dexseq, "gene_full_by_tp", tp)
    gene_sig <- get_tp_table(dexseq, "gene_sig_by_tp", tp)
    add_row(summary_row(
      method = "DEXSeq",
      tp = tp,
      features = "exon bins",
      significance_metric = "BH-adjusted P value",
      significance_cutoff = dexseq$meta$padj_cutoff %||% 0.10,
      n_features_tested = count_tested_rows(exon_full, c("pvalue", "padj")),
      n_sig_features = nrow(exon_sig),
      n_genes_tested = if (nrow(gene_full) > 0) {
        count_tested_rows(gene_full, c("pvalue", "padj"))
      } else {
        count_unique_genes(
          exon_full,
          tested_only = TRUE,
          statistic_candidates = c("pvalue", "padj")
        )
      },
      n_sig_genes = nrow(gene_sig),
      notes = "Significant genes contain at least one significant exon-bin usage test."
    ))

    # DRIMSeq: transcript features are tested within genes and confirmed with stageR.
    tx_full <- get_tp_table(drimseq, "tx_full_by_tp", tp)
    tx_sig <- get_tp_table(drimseq, "tx_sig_by_tp", tp)
    gene_full <- get_tp_table(drimseq, "gene_full_by_tp", tp)
    gene_sig <- get_tp_table(drimseq, "gene_sig_by_tp", tp)
    add_row(summary_row(
      method = "DRIMSeq",
      tp = tp,
      features = "transcripts",
      significance_metric = "stageR adjusted P value",
      significance_cutoff = drimseq$meta$padj_cutoff %||% 0.10,
      n_features_tested = count_tested_rows(tx_full, c("pvalue", "adj_pvalue", "padj")),
      n_sig_features = nrow(tx_sig),
      n_genes_tested = if (nrow(gene_full) > 0) {
        count_tested_rows(gene_full, c("pvalue", "adj_pvalue", "padj"))
      } else {
        count_unique_genes(
          tx_full,
          tested_only = TRUE,
          statistic_candidates = c("pvalue", "adj_pvalue", "padj")
        )
      },
      n_sig_genes = nrow(gene_sig),
      notes = "Transcript-level confirmation is conditional on gene-level screening."
    ))

    # SUPPA2 diffSplice reports nominal P values; do not mislabel this as padj.
    event_full <- get_tp_table(suppa, "event_full_by_tp", tp)
    event_sig <- get_tp_table(suppa, "event_sig_by_tp", tp)
    gene_sig <- get_tp_table(suppa, "gene_sig_by_tp", tp)
    add_row(summary_row(
      method = "SUPPA2",
      tp = tp,
      features = "splicing events",
      significance_metric = "unadjusted diffSplice P value",
      significance_cutoff = suppa$meta$p_cutoff %||% 0.10,
      n_features_tested = count_tested_rows(event_full, c("pvalue")),
      n_sig_features = nrow(event_sig),
      n_genes_tested = count_unique_genes(
        event_full,
        tested_only = TRUE,
        statistic_candidates = c("pvalue")
      ),
      n_sig_genes = nrow(gene_sig),
      notes = "SUPPA2 diffSplice does not provide an adjusted-P column in this pipeline; the manuscript threshold is nominal P < 0.10."
    ))

    # IsoformSwitchAnalyzeR: use the same strict cache and call definition as Figure 6.
    switches <- as.data.frame(
      fig6_analysis$tables$isoform_switches %||% data.frame(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    tp_col <- first_present_column(switches, c("timepoint", "condition_2"))
    sig_tp <- if (!is.null(tp_col)) {
      switches[as.character(switches[[tp_col]]) == tp, , drop = FALSE]
    } else {
      data.frame()
    }

    background <- as.data.frame(
      fig6_analysis$tables$background_transcripts %||% strict_switch_list$isoformFeatures %||% data.frame(),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    bg_tp_col <- first_present_column(background, c("timepoint", "condition_2"))
    if (!is.null(bg_tp_col)) {
      background <- background[as.character(background[[bg_tp_col]]) == tp, , drop = FALSE]
    }
    tx_col <- first_present_column(background, c(
      "isoform_id", "transcript_id", "feature_id",
      "gencode_transcript_id"
    ))
    gene_col <- first_present_column(background, c(
      "gene_id", "ensgene", "gencode_gene_id"
    ))
    n_bg_tx <- if (is.null(tx_col)) NA_integer_ else {
      as.integer(length(unique(background[[tx_col]][
        !is.na(background[[tx_col]]) & nzchar(as.character(background[[tx_col]]))
      ])))
    }
    n_bg_gene <- if (is.null(gene_col)) NA_integer_ else {
      as.integer(length(unique(strip_ens_version(background[[gene_col]][
        !is.na(background[[gene_col]]) & nzchar(as.character(background[[gene_col]]))
      ]))))
    }
    sig_tx_col <- first_present_column(sig_tp, c(
      "isoform_id", "transcript_id", "feature_id",
      "gencode_transcript_id"
    ))
    sig_gene_col <- first_present_column(sig_tp, c(
      "gene_id", "ensgene", "gencode_gene_id"
    ))
    n_sig_tx <- if (is.null(sig_tx_col)) 0L else {
      as.integer(length(unique(sig_tp[[sig_tx_col]][
        !is.na(sig_tp[[sig_tx_col]]) & nzchar(as.character(sig_tp[[sig_tx_col]]))
      ])))
    }
    n_sig_gene <- if (is.null(sig_gene_col)) 0L else {
      as.integer(length(unique(strip_ens_version(sig_tp[[sig_gene_col]][
        !is.na(sig_tp[[sig_gene_col]]) & nzchar(as.character(sig_tp[[sig_gene_col]]))
      ]))))
    }
    add_row(summary_row(
      method = "IsoformSwitchAnalyzeR",
      tp = tp,
      features = "transcripts",
      significance_metric = "isoform-switch q value",
      significance_cutoff = 0.05,
      effect_size_cutoff = "|dIF| > 0.10",
      n_features_tested = n_bg_tx,
      n_sig_features = n_sig_tx,
      n_genes_tested = n_bg_gene,
      n_sig_genes = n_sig_gene,
      notes = "Counts use the strict Figure 6 cache and the same q-value/dIF call definition."
    ))
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# ------------------------------------------------------------------------------
# Per-sample sequencing and mapping QC
# ------------------------------------------------------------------------------
read_sample_metadata <- function() {
  path <- first_existing(
    c(
      Sys.getenv("SAMPLE_METADATA_TSV"),
      file.path(PROJECT_ROOT, "config", "metadata", "samples.tsv")
    ),
    required = FALSE,
    label = "sample metadata TSV"
  )
  if (is.na(path)) return(data.frame())

  x <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  sample_col <- first_present_column(x, c("sample_id", "names", "sample", "Sample"))
  if (is.null(sample_col)) return(data.frame())
  x$sample_id <- as.character(x[[sample_col]])
  if ("exclude" %in% colnames(x)) {
    excluded <- toupper(as.character(x$exclude)) %in% c("TRUE", "T", "1", "YES")
    x <- x[!excluded, , drop = FALSE]
  }
  x$sample_key <- normalize_sample_key(x$sample_id)
  x
}

read_salmon_meta <- function(path) {
  if (!file.exists(path)) return(list())
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    warning(
      "Package 'jsonlite' is not installed; Salmon meta_info.json metrics ",
      "will be omitted. Install jsonlite or use a MultiQC export."
    )
    return(list())
  }
  jsonlite::fromJSON(path, simplifyVector = TRUE)
}

collect_salmon_qc <- function() {
  coldata_path <- first_existing(
    c(
      Sys.getenv("SALMON_COLDATA_TSV"),
      file.path(PROJECT_ROOT, "results", "salmon", "coldata.tsv"),
      file.path(PROJECT_ROOT, "results", "counts", "salmon", "coldata.tsv")
    ),
    required = FALSE,
    label = "Salmon coldata.tsv"
  )

  records <- list()
  if (!is.na(coldata_path)) {
    cd <- utils::read.delim(coldata_path, stringsAsFactors = FALSE, check.names = FALSE)
    sample_col <- first_present_column(cd, c("names", "sample_id", "sample"))
    file_col <- first_present_column(cd, c("files", "quant_file", "quant.sf"))
    if (!is.null(sample_col) && !is.null(file_col)) {
      for (i in seq_len(nrow(cd))) {
        records[[length(records) + 1L]] <- data.frame(
          sample_id = as.character(cd[[sample_col]][[i]]),
          condition_salmon = if ("condition" %in% colnames(cd)) as.character(cd$condition[[i]]) else NA_character_,
          quant_file = resolve_relative_path(as.character(cd[[file_col]][[i]])),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(records) == 0) {
    salmon_root <- first_existing(
      c(
        Sys.getenv("SALMON_QC_DIR"),
        file.path(PROJECT_ROOT, "results", "salmon"),
        file.path(PROJECT_ROOT, "results", "counts", "salmon")
      ),
      required = FALSE,
      label = "Salmon QC directory"
    )
    if (!is.na(salmon_root) && dir.exists(salmon_root)) {
      quant_files <- list.files(
        salmon_root,
        pattern = "^quant\\.sf$",
        recursive = TRUE,
        full.names = TRUE
      )
      records <- lapply(quant_files, function(path) {
        data.frame(
          sample_id = basename(dirname(path)),
          condition_salmon = basename(dirname(dirname(path))),
          quant_file = normalizePath(path, mustWork = FALSE),
          stringsAsFactors = FALSE
        )
      })
    }
  }

  if (length(records) == 0) return(data.frame())
  records <- do.call(rbind, records)

  out <- lapply(seq_len(nrow(records)), function(i) {
    quant_file <- records$quant_file[[i]]
    quant <- if (file.exists(quant_file)) {
      utils::read.delim(quant_file, stringsAsFactors = FALSE, check.names = FALSE)
    } else {
      data.frame()
    }
    meta_path <- file.path(dirname(quant_file), "aux_info", "meta_info.json")
    meta <- read_salmon_meta(meta_path)

    data.frame(
      sample_id = records$sample_id[[i]],
      sample_key = normalize_sample_key(records$sample_id[[i]]),
      condition_salmon = records$condition_salmon[[i]],
      salmon_processed_fragments = as_numeric_or_na(meta$num_processed %||% NA),
      salmon_mapped_fragments = as_numeric_or_na(meta$num_mapped %||% NA),
      salmon_mapping_pct = as_numeric_or_na(meta$percent_mapped %||% NA),
      salmon_quantified_fragments = if (
        nrow(quant) > 0 && "NumReads" %in% colnames(quant)
      ) {
        sum(suppressWarnings(as.numeric(quant$NumReads)), na.rm = TRUE)
      } else {
        NA_real_
      },
      salmon_targets_in_quant = nrow(quant),
      salmon_version = as.character(meta$salmon_version %||% NA_character_),
      salmon_meta_found = file.exists(meta_path),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

star_value <- function(parsed, key) {
  hit <- parsed$value[parsed$key == key]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

collect_star_qc <- function() {
  roots <- unique(c(
    Sys.getenv("STAR_QC_DIR"),
    file.path(PROJECT_ROOT, "results", "star"),
    file.path(PROJECT_ROOT, "results", "counts", "star"),
    file.path(PROJECT_ROOT, "results", "alignment", "star"),
    file.path(PROJECT_ROOT, "results", "alignments", "star")
  ))
  roots <- roots[nzchar(roots) & dir.exists(roots)]

  logs <- unique(unlist(lapply(roots, function(root) {
    list.files(
      root,
      pattern = "Log\\.final\\.out$",
      recursive = TRUE,
      full.names = TRUE
    )
  }), use.names = FALSE))
  if (length(logs) == 0) return(data.frame())

  rows <- lapply(logs, function(path) {
    lines <- readLines(path, warn = FALSE)
    parts <- strsplit(lines[grepl("\\|", lines)], "\\|")
    parsed <- data.frame(
      key = trimws(vapply(parts, `[[`, character(1), 1)),
      value = trimws(vapply(parts, function(z) paste(z[-1], collapse = "|"), character(1))),
      stringsAsFactors = FALSE
    )

    filename <- basename(path)
    sample_id <- if (identical(filename, "Log.final.out")) {
      basename(dirname(path))
    } else {
      sub("[._-]?Log\\.final\\.out$", "", filename)
    }

    unique_pct <- as_numeric_or_na(star_value(parsed, "Uniquely mapped reads %"))
    multi_pct <- as_numeric_or_na(star_value(parsed, "% of reads mapped to multiple loci"))
    too_many_pct <- as_numeric_or_na(star_value(parsed, "% of reads mapped to too many loci"))

    mapped_components <- c(unique_pct, multi_pct, too_many_pct)
    total_mapped_pct <- if (all(is.na(mapped_components))) {
      NA_real_
    } else {
      sum(mapped_components, na.rm = TRUE)
    }

    data.frame(
      sample_id_star = sample_id,
      sample_key = normalize_sample_key(sample_id),
      star_input_fragments = as_numeric_or_na(star_value(parsed, "Number of input reads")),
      star_uniquely_mapped_fragments = as_numeric_or_na(star_value(parsed, "Uniquely mapped reads number")),
      star_uniquely_mapped_pct = unique_pct,
      star_multimapped_fragments = as_numeric_or_na(star_value(parsed, "Number of reads mapped to multiple loci")),
      star_multimapped_pct = multi_pct,
      star_too_many_loci_pct = too_many_pct,
      star_total_mapped_pct = total_mapped_pct,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  if (anyDuplicated(out$sample_key)) {
    warning(
      "Multiple STAR Log.final.out files resolved to the same sample ID; ",
      "using the first file for each sample."
    )
    out <- out[!duplicated(out$sample_key), , drop = FALSE]
  }
  out
}

build_sample_qc_table <- function() {
  metadata <- read_sample_metadata()
  salmon_qc <- collect_salmon_qc()
  star_qc <- collect_star_qc()

  keys <- unique(c(
    metadata$sample_key %||% character(),
    salmon_qc$sample_key %||% character(),
    star_qc$sample_key %||% character()
  ))
  keys <- keys[!is.na(keys) & nzchar(keys)]
  if (length(keys) == 0) {
    return(data.frame(
      note = paste(
        "No Salmon or STAR QC files were found.",
        "Provide Salmon quant directories containing aux_info/meta_info.json",
        "and STAR Log.final.out files, or set SALMON_QC_DIR/STAR_QC_DIR."
      )
    ))
  }

  rows <- lapply(keys, function(key) {
    m <- if (nrow(metadata) > 0) metadata[metadata$sample_key == key, , drop = FALSE] else data.frame()
    s <- if (nrow(salmon_qc) > 0) salmon_qc[salmon_qc$sample_key == key, , drop = FALSE] else data.frame()
    z <- if (nrow(star_qc) > 0) star_qc[star_qc$sample_key == key, , drop = FALSE] else data.frame()

    sample_id <- if (nrow(m) > 0) {
      m$sample_id[[1]]
    } else if (nrow(s) > 0) {
      s$sample_id[[1]]
    } else {
      z$sample_id_star[[1]]
    }
    condition_col <- if (nrow(m) > 0) {
      first_present_column(m, c("condition", "timepoint", "group"))
    } else {
      NULL
    }
    condition <- if (!is.null(condition_col)) {
      as.character(m[[condition_col]][[1]])
    } else if (nrow(s) > 0) {
      as.character(s$condition_salmon[[1]])
    } else {
      NA_character_
    }
    layout_col <- if (nrow(m) > 0) {
      first_present_column(
        m,
        c("layout", "library_layout", "sequencing_layout")
      )
    } else {
      NULL
    }
    layout <- if (!is.null(layout_col)) {
      as.character(m[[layout_col]][[1]])
    } else {
      NA_character_
    }
    replicate_col <- if (nrow(m) > 0) {
      first_present_column(m, c("replicate", "rep", "Replicate"))
    } else {
      NULL
    }
    replicate <- if (!is.null(replicate_col)) {
      as.character(m[[replicate_col]][[1]])
    } else {
      NA_character_
    }
    if ((is.na(replicate) || !nzchar(replicate)) && grepl("[_-][0-9]+$", sample_id)) {
      replicate <- sub("^.*[_-]([0-9]+)$", "\\1", sample_id)
    }

    star_input <- if (nrow(z) > 0) z$star_input_fragments[[1]] else NA_real_
    salmon_input <- if (nrow(s) > 0) s$salmon_processed_fragments[[1]] else NA_real_
    salmon_quantified <- if (nrow(s) > 0) {
      s$salmon_quantified_fragments[[1]]
    } else {
      NA_real_
    }
    reported_input_count <- if (is.finite(salmon_input)) {
      salmon_input
    } else if (is.finite(salmon_quantified)) {
      salmon_quantified
    } else {
      star_input
    }
    count_definition <- if (is.finite(salmon_input)) {
      "Salmon processed fragments/read pairs"
    } else if (is.finite(salmon_quantified)) {
      "Sum of Salmon estimated NumReads (quantified fragments; not total input)"
    } else if (is.finite(star_input)) {
      "STAR Number of input reads for the pooled alignment unit"
    } else {
      NA_character_
    }

    condition_codes <- c("C1", "H1", "H3", "H24")
    condition_from_star <- if (
      nrow(z) > 0 && toupper(as.character(z$sample_id_star[[1]])) %in% condition_codes
    ) {
      toupper(as.character(z$sample_id_star[[1]]))
    } else {
      NA_character_
    }
    if (is.na(condition) || !nzchar(condition)) condition <- condition_from_star

    qc_level <- if (
      nrow(z) > 0 && nrow(s) == 0 && nrow(m) == 0 &&
        !is.na(condition_from_star)
    ) {
      "condition-level pooled alignment"
    } else {
      "sample"
    }
    metric_sources <- paste(
      c(if (nrow(s) > 0) "Salmon", if (nrow(z) > 0) "STAR"),
      collapse = "; "
    )

    notes <- character()
    if (nrow(s) == 0) notes <- c(notes, "Salmon QC unavailable")
    if (nrow(s) > 0 && !isTRUE(s$salmon_meta_found[[1]])) {
      notes <- c(notes, "Salmon meta_info.json missing; quant.sf alone cannot provide mapping rate")
    }
    if (nrow(z) == 0) notes <- c(notes, "STAR Log.final.out unavailable")
    if (
      is.finite(star_input) && is.finite(salmon_input) &&
        abs(star_input - salmon_input) / max(star_input, salmon_input) > 0.001
    ) {
      notes <- c(notes, "STAR and Salmon processed-fragment totals differ")
    }
    if (identical(qc_level, "condition-level pooled alignment")) {
      notes <- c(
        notes,
        "STAR metrics apply to the pooled condition alignment, not to individual biological replicates"
      )
    }

    data.frame(
      qc_unit_id = sample_id,
      qc_level = qc_level,
      condition = condition,
      replicate = replicate,
      library_layout = layout,
      metric_sources = metric_sources,
      reported_input_count = reported_input_count,
      reported_input_count_definition = count_definition,
      star_uniquely_mapped_fragments = if (nrow(z) > 0) z$star_uniquely_mapped_fragments[[1]] else NA_real_,
      star_uniquely_mapped_pct = if (nrow(z) > 0) z$star_uniquely_mapped_pct[[1]] else NA_real_,
      star_multimapped_pct = if (nrow(z) > 0) z$star_multimapped_pct[[1]] else NA_real_,
      star_total_mapped_pct = if (nrow(z) > 0) z$star_total_mapped_pct[[1]] else NA_real_,
      salmon_processed_fragments = salmon_input,
      salmon_mapped_fragments = if (nrow(s) > 0) s$salmon_mapped_fragments[[1]] else NA_real_,
      salmon_mapping_pct = if (nrow(s) > 0) s$salmon_mapping_pct[[1]] else NA_real_,
      salmon_quantified_fragments = salmon_quantified,
      salmon_targets_in_quant = if (nrow(s) > 0) s$salmon_targets_in_quant[[1]] else NA_integer_,
      salmon_version = if (nrow(s) > 0) s$salmon_version[[1]] else NA_character_,
      qc_notes = if (length(notes) == 0) NA_character_ else paste(notes, collapse = "; "),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  tp_order <- match(out$condition, c("C1", "H1", "H3", "H24"))
  level_order <- match(out$qc_level, c("sample", "condition-level pooled alignment"))
  out <- out[
    order(tp_order, level_order, out$replicate, out$qc_unit_id, na.last = TRUE),
    ,
    drop = FALSE
  ]
  rownames(out) <- NULL
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
    } else NULL,

    # Fallback if the standalone RDS stores complete_tbl at its top level.
    fig1_go_top_level_complete = if (!is.null(fig1_go)) {
      fig1_go$complete_tbl
    } else NULL,

    # Figure 1 gene-level object produced by the updated plotting script.
    fig1_pathway_enrichment_complete = if (!is.null(fig1)) {
      fig1$pathway_enrichment$enrichment$complete_tbl
    } else NULL,

    # Fallback for an alternate nested structure.
    fig1_enrichment_complete = if (!is.null(fig1)) {
      fig1$enrichment$complete_tbl
    } else NULL
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
    if (length(genes) == 0) return(NULL)
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

if (!is.null(fig5_enrichment)) {
  if (!is.null(fig5_enrichment$enrichment$as_pathways_tbl) && nrow(fig5_enrichment$enrichment$as_pathways_tbl) > 0) {
    pathway_asg <- fig5_enrichment$enrichment$as_pathways_tbl
  }
  if (!is.null(fig5_enrichment$enrichment$deg_pathways_tbl) && nrow(fig5_enrichment$enrichment$deg_pathways_tbl) > 0) {
    pathway_deg <- fig5_enrichment$enrichment$deg_pathways_tbl
  }
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

method_summary <- make_standard_method_summary(
  deseq = deseq,
  dexseq = dexseq,
  drimseq = drimseq,
  suppa = suppa,
  fig6_analysis = fig6_analysis,
  strict_switch_list = strict_switch_list,
  timepoints = timepoints,
  reference_condition = "C1"
)

message_step("Collecting sample- and condition-level sequencing QC")
sample_qc <- build_sample_qc_table()

sheet_titles <- c(
  "01 README & Contents",
  "02 Methods Summary",
  "03 DESeq2 Genes",
  "04 DEXSeq Genes",
  "05 DEXSeq Exon Bins",
  "06 DRIMSeq Genes",
  "07 DRIMSeq Transcripts",
  "08 SUPPA2 Genes",
  "09 SUPPA2 Events",
  "10 Integrated Status (Fig 4B)",
  "11 Transitions (Fig 4C)",
  "12 DEG Clusters (Fig 1)",
  "13 Cluster GO (Fig 1C)",
  "14 Pathway Summary (Fig 5A)",
  "15 ASG Pathway Enrichment",
  "16 DEG Pathway Enrichment",
  "17 Isoform Switches (Fig 6)",
  "18 Isoform Biotypes",
  "19 DRIMSeq Biotypes",
  "20 NMD Switch Pairs",
  "21 NMD Statistics",
  "22 PTC Enrichment",
  "23 Sequencing QC"
)

sheet_descriptions <- c(
  "Workbook overview, analysis definitions, thresholds, provenance notes, and a table of contents describing every worksheet.",
  "Standardized method-by-timepoint summary of the feature type, significance rule, number of features and genes tested, and number called significant.",
  "Significant gene-level differential-expression results from DESeq2 for all hypoxia-versus-C1 contrasts.",
  "DEXSeq gene-level summaries for genes containing at least one significant exon-bin usage change.",
  "Significant DEXSeq exon-bin usage results for all hypoxia-versus-C1 contrasts.",
  "Significant DRIMSeq gene-level differential-transcript-usage results.",
  "Significant DRIMSeq transcript-level results annotated with GENCODE v45 transcript and gene biotypes.",
  "Unique genes containing at least one significant SUPPA2 splicing event.",
  "Significant SUPPA2 event-level results, including event class, delta PSI, nominal P value, and direction.",
  "Union of genes detected by DESeq2, DEXSeq, DRIMSeq, or SUPPA2, with method- and timepoint-specific status and the classifications plotted in Figure 4B.",
  "Genes belonging to the highlighted temporal regulatory-transition groups plotted in Figure 4C.",
  "Figure 1 temporal DEG-cluster membership together with scaled mean expression values used to define the clusters.",
  "Complete significant cluster-level GO enrichment results underlying the representative terms displayed in Figure 1C.",
  "Counts of differential-expression and alternative-splicing pathway representation underlying Figure 5A.",
  "Complete pathway-enrichment results for alternatively spliced genes.",
  "Complete pathway-enrichment results for differentially expressed genes.",
  "Significant IsoformSwitchAnalyzeR calls from the strict analysis cache underlying Figure 6, with transcript biotype and NMD/PTC annotations where available.",
  "Transcript-biotype enrichment among IsoformSwitchAnalyzeR calls relative to all eligible expressed transcripts in each contrast.",
  "Transcript-biotype enrichment among significant DRIMSeq transcripts relative to all transcripts tested in each contrast.",
  "Paired gained/lost isoforms used to evaluate NMD-sensitive switch consequences at each hypoxic timepoint.",
  "IsoformSwitchAnalyzeR NMD consequence counts and enrichment statistics.",
  "Fisher exact tests of PTC-positive transcripts among significant switch transcripts relative to the PTC-evaluable expressed-transcript background.",
  "Sequencing and mapping QC. Salmon metrics are reported per biological sample; STAR uniquely mapped-read metrics are reported at the condition-pooled level when that was the alignment unit."
)

if (
  length(sheet_titles) != length(sheet_descriptions) ||
    anyDuplicated(sheet_titles) ||
    any(nchar(sheet_titles) > 31)
) {
  stop(
    "Worksheet titles/descriptions are inconsistent, duplicated, or exceed ",
    "Excel's 31-character worksheet-name limit."
  )
}

general_information <- data.frame(
  section = "Workbook information",
  item = c(
    "Manuscript title",
    "Workbook purpose",
    "Contrasts",
    "Significance thresholds",
    "Feature-count definitions",
    "Annotation",
    "Biotype enrichment background",
    "NMD consequence analysis",
    "PTC background enrichment",
    "Sequencing QC interpretation",
    "Strict IsoformSwitchAnalyzeR source",
    "Data availability",
    "Generated on"
  ),
  description = c(
    "Gene-level RNA-seq obscures extensive transcript-level remodeling during hypoxic adaptation",
    "Curated supplementary results supporting differential expression, alternative splicing and transcript usage, isoform switching, transcript-biotype enrichment, integrated gene classification, temporal transitions, pathway analyses, and sequencing QC.",
    paste(paste(names(timepoint_labels), timepoint_labels, sep = " = "), collapse = "; "),
    "DESeq2 and DEXSeq use BH-adjusted P < 0.10; DRIMSeq uses stageR-adjusted P < 0.10; SUPPA2 uses the nominal diffSplice P < 0.10 produced by this pipeline; IsoformSwitchAnalyzeR uses q < 0.05 and |dIF| > 0.10. Tool-specific rules are tabulated in worksheet 02.",
    "A tested feature is a gene for DESeq2, an exon bin for DEXSeq, a transcript for DRIMSeq and IsoformSwitchAnalyzeR, and a splicing event for SUPPA2. The number of tested genes is also reported so feature- and gene-level totals are not conflated.",
    "Transcript and gene biotypes are taken from the same GENCODE v45 basic GTF used for project quantification.",
    "For each C1-to-hypoxia contrast, the IsoformSwitchAnalyzeR background is all eligible unique transcripts in isoformFeatures and the DRIMSeq background is all unique transcripts in tx_full_by_tp. GENCODE supplies annotation labels but does not define the tested background.",
    "NMD status was recalculated from the strict cached switch calls and their PTC annotations using analyzeSwitchConsequences(consequencesToAnalyze = 'NMD_status', alpha = 0.05, dIFcutoff = 0.10, onlySigIsoforms = FALSE). An NMD-sensitive gain denotes a PTC-positive isoform used more in hypoxia paired with a PTC-negative isoform used less.",
    "For each timepoint and switch direction, PTC-positive switching transcripts were compared with other PTC-evaluable transcripts in the strict isoformFeatures background using two-sided Fisher exact tests, with BH correction across timepoint-by-direction tests.",
    "Salmon num_processed and num_mapped are fragment/read-pair counts from each replicate's aux_info/meta_info.json; sum(NumReads) from quant.sf is retained as the quantified-fragment total. STAR 'Number of input reads' and uniquely mapped fractions are condition-level because the project STAR alignment was performed after pooling by condition. These pooled values are not assigned to individual replicates.",
    basename(strict_switch_rds),
    "The workbook contains significant method-level calls and complete enrichment outputs used for the reported analyses. Full unfiltered statistical results remain in the cached RDS analysis objects.",
    as.character(Sys.time())
  ),
  stringsAsFactors = FALSE
)

worksheet_guide <- data.frame(
  section = "Worksheet guide",
  item = sheet_titles,
  description = sheet_descriptions,
  stringsAsFactors = FALSE
)
readme <- rbind(general_information, worksheet_guide)

sheet_tables <- stats::setNames(
  list(
    readme,
    method_summary,
    deseq_degs,
    dexseq_genes,
    dexseq_exons,
    drimseq_genes,
    drimseq_tx,
    suppa_genes,
    suppa_events,
    integrated,
    transitions,
    fig1_clusters,
    fig1_go_full,
    pathway_representation,
    pathway_asg,
    pathway_deg,
    isoform_switches,
    isoform_biotype_enrichment,
    drimseq_biotype_enrichment,
    isoform_nmd_details,
    isoform_nmd_stats,
    isoform_ptc_summary,
    sample_qc
  ),
  sheet_titles
)

for (sheet_name in names(sheet_tables)) {
  add_sheet(wb, sheet_name, sheet_tables[[sheet_name]])
}

openxlsx::saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
message_step("Done: ", OUT_XLSX)
