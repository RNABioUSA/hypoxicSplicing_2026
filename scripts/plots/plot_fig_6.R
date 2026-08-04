# ============================================================
# plot_fig_6.R
# ============================================================
#
# Refactored Figure 6 workflow for IsoformSwitchAnalyzeR.
#
# The script assumes that the expensive sequence/consequence analyses have
# already been completed and cached in switchList_filtered_analyzed.rds. It does
# not rerun isoformSwitchTestDEXSeq(), CPC2, PFAM, SignalP, IUPred2A, DeepTMHMM,
# DeepLoc2, or analyzeSwitchConsequences().
#
# Public functions:
#   get_gencode_biotype_annotation()
#   annotate_transcript_biotypes()
#   annotate_gene_biotypes()
#   find_isoform_switch_cache()
#   run_isoform_switch_analysis()
#   plot_switch_overlap_and_deg_fraction()       # Figure 6A
#   plot_isoform_switch_examples()                # Figure 6B
#   plot_switch_usage_vs_gene_expression()        # Figure 6C
#   plot_switch_consequence_summary()             # Figure 6D
#   plot_fig_6_all()
#
# Example:
#   source(file.path(PROJECT_ROOT, "scripts/plots/plot_fig_6.R"))
#   fig6 <- plot_fig_6_all()
#

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
if (!nzchar(PROJECT_ROOT)) stop("PROJECT_ROOT env var not set.")

HELPERS_FILE <- file.path(PROJECT_ROOT, "scripts", "utils", "helpers.R")
if (file.exists(HELPERS_FILE)) source(HELPERS_FILE)

COLORS_FILE <- file.path(PROJECT_ROOT, "scripts", "utils", "color_palette.R")
if (file.exists(COLORS_FILE)) source(COLORS_FILE)

`%||%` <- function(x, y) if (is.null(x)) y else x

.fig6_analysis_version <- "2026-07-28.1"

.fig6_check_pkg <- function(pkg) {
  if (exists(".check_pkg", mode = "function")) {
    return(.check_pkg(pkg))
  }
  for (p in pkg) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Package '", p, "' is required.")
    }
  }
  invisible(TRUE)
}

.fig6_strip_version <- function(x) {
  if (exists(".strip_ens_version", mode = "function")) {
    return(.strip_ens_version(x))
  }
  sub("\\.[0-9]+$", "", as.character(x))
}

.fig6_first_existing <- function(paths, required = TRUE, label = "file") {
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  hit <- paths[file.exists(paths)]
  if (length(hit) > 0) {
    return(normalizePath(hit[[1]], mustWork = TRUE))
  }
  if (isTRUE(required)) {
    stop("Could not locate ", label, ". Tried:\n  - ", paste(paths, collapse = "\n  - "))
  }
  NA_character_
}

.fig6_first_col <- function(df, candidates, required = FALSE, label = "column") {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) {
    return(hit[[1]])
  }
  if (isTRUE(required)) {
    stop(
      "Could not identify ", label, ". Expected one of: ",
      paste(candidates, collapse = ", "),
      ". Available columns: ", paste(colnames(df), collapse = ", ")
    )
  }
  NULL
}

.fig6_bind_rows <- function(dfs) {
  dfs <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, dfs)
  if (length(dfs) == 0) {
    return(data.frame())
  }
  cols <- unique(unlist(lapply(dfs, names), use.names = FALSE))
  dfs <- lapply(dfs, function(x) {
    for (nm in setdiff(cols, names(x))) x[[nm]] <- NA
    x[, cols, drop = FALSE]
  })
  out <- do.call(rbind, dfs)
  rownames(out) <- NULL
  out
}

.fig6_make_unique <- function(df, key) {
  if (is.null(df) || nrow(df) == 0) {
    return(df)
  }
  df[!duplicated(df[[key]]) & !is.na(df[[key]]) & df[[key]] != "", , drop = FALSE]
}

.fig6_outdir <- function(outdir = NULL) {
  if (!is.null(outdir) && nzchar(outdir)) {
    resolved <- outdir
  } else if (exists(".get_results_dir", mode = "function")) {
    resolved <- file.path(.get_results_dir(), "plots", "fig_6")
  } else {
    resolved <- file.path(PROJECT_ROOT, "results", "plots", "fig_6")
  }
  dir.create(resolved, recursive = TRUE, showWarnings = FALSE)
  normalizePath(resolved, mustWork = FALSE)
}

# --------------------------------------------------
# GENCODE transcript and gene biotypes
# --------------------------------------------------

.fig6_find_gtf <- function(
  gtf_path = NULL,
  env_file = file.path(PROJECT_ROOT, "config", "environment", "analysis_environment.env"),
  resources_dir = file.path(PROJECT_ROOT, "resources")
) {
  if (!is.null(gtf_path) && nzchar(gtf_path)) {
    return(.fig6_first_existing(gtf_path, label = "GENCODE GTF"))
  }

  if (
    exists(".resolve_ref_from_env", mode = "function") &&
      file.exists(env_file) &&
      dir.exists(resources_dir)
  ) {
    ref <- tryCatch(
      .resolve_ref_from_env(env_file = env_file, resources_dir = resources_dir),
      error = function(e) NULL
    )
    if (!is.null(ref) && file.exists(ref$gtf_path)) {
      return(normalizePath(ref$gtf_path, mustWork = TRUE))
    }
  }

  candidates <- if (dir.exists(resources_dir)) {
    list.files(
      resources_dir,
      pattern = "\\.gtf(\\.gz)?$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
  } else {
    character()
  }

  # Prefer GENCODE v45 and then "basic" annotation if multiple files exist.
  if (length(candidates) > 1) {
    score <- 4L * grepl("gencode", basename(candidates), ignore.case = TRUE) +
      2L * grepl("v45|release[_-]?45", basename(candidates), ignore.case = TRUE) +
      1L * grepl("basic", basename(candidates), ignore.case = TRUE)
    candidates <- candidates[order(score, decreasing = TRUE)]
  }

  .fig6_first_existing(candidates, label = "GENCODE GTF under PROJECT_ROOT/resources")
}

get_gencode_biotype_annotation <- function(
  gtf_path = NULL,
  env_file = file.path(PROJECT_ROOT, "config", "environment", "analysis_environment.env"),
  resources_dir = file.path(PROJECT_ROOT, "resources"),
  cache_file = file.path(resources_dir, "cache", "gencode_transcript_biotypes.rds"),
  force = FALSE
) {
  .fig6_check_pkg(c("rtracklayer", "S4Vectors"))
  gtf_path <- .fig6_find_gtf(
    gtf_path = gtf_path,
    env_file = env_file,
    resources_dir = resources_dir
  )

  gtf_info <- file.info(gtf_path)
  source_stamp <- list(
    path = normalizePath(gtf_path, mustWork = TRUE),
    size = unname(gtf_info$size),
    mtime = as.character(gtf_info$mtime)
  )

  if (!isTRUE(force) && file.exists(cache_file)) {
    cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (
      is.list(cached) &&
        identical(cached$meta$source, source_stamp) &&
        is.data.frame(cached$transcripts)
    ) {
      return(cached)
    }
  }

  message("[Figure 6] Importing transcript biotypes from: ", gtf_path)
  gtf <- rtracklayer::import(gtf_path)
  meta <- as.data.frame(S4Vectors::mcols(gtf), stringsAsFactors = FALSE)

  if (!"type" %in% colnames(meta)) {
    stop("The imported GTF does not contain a feature 'type' column.")
  }
  tx <- meta[as.character(meta$type) == "transcript", , drop = FALSE]
  if (nrow(tx) == 0) stop("No transcript features were found in the GTF.")

  transcript_id_col <- .fig6_first_col(
    tx,
    c("transcript_id", "transcriptId"),
    required = TRUE,
    label = "GENCODE transcript ID"
  )
  gene_id_col <- .fig6_first_col(
    tx,
    c("gene_id", "geneId"),
    required = TRUE,
    label = "GENCODE gene ID"
  )
  tx_type_col <- .fig6_first_col(
    tx,
    c("transcript_type", "transcript_biotype"),
    required = TRUE,
    label = "GENCODE transcript biotype"
  )
  gene_type_col <- .fig6_first_col(
    tx,
    c("gene_type", "gene_biotype"),
    required = FALSE
  )
  gene_name_col <- .fig6_first_col(tx, c("gene_name", "gene"), required = FALSE)
  tx_name_col <- .fig6_first_col(tx, c("transcript_name"), required = FALSE)
  tsl_col <- .fig6_first_col(tx, c("transcript_support_level"), required = FALSE)

  transcript_annot <- data.frame(
    transcript_id_full = as.character(tx[[transcript_id_col]]),
    transcript_id = .fig6_strip_version(tx[[transcript_id_col]]),
    gene_id_full = as.character(tx[[gene_id_col]]),
    gene_id = .fig6_strip_version(tx[[gene_id_col]]),
    gene_symbol = if (is.null(gene_name_col)) NA_character_ else as.character(tx[[gene_name_col]]),
    transcript_name = if (is.null(tx_name_col)) NA_character_ else as.character(tx[[tx_name_col]]),
    transcript_biotype = as.character(tx[[tx_type_col]]),
    gene_biotype = if (is.null(gene_type_col)) NA_character_ else as.character(tx[[gene_type_col]]),
    transcript_support_level = if (is.null(tsl_col)) NA_character_ else as.character(tx[[tsl_col]]),
    stringsAsFactors = FALSE
  )
  transcript_annot <- .fig6_make_unique(transcript_annot, "transcript_id")

  gene_annot <- transcript_annot[
    !duplicated(transcript_annot$gene_id),
    c("gene_id_full", "gene_id", "gene_symbol", "gene_biotype"),
    drop = FALSE
  ]

  out <- list(
    meta = list(
      source = source_stamp,
      n_transcripts = nrow(transcript_annot),
      n_genes = nrow(gene_annot)
    ),
    transcripts = transcript_annot,
    genes = gene_annot
  )

  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(out, cache_file)
  out
}

annotate_transcript_biotypes <- function(df, annotation, transcript_col = NULL) {
  if (is.null(df) || nrow(df) == 0) {
    return(df)
  }
  if (is.null(annotation$transcripts)) stop("Annotation object lacks $transcripts.")

  transcript_col <- transcript_col %||% .fig6_first_col(
    df,
    c("transcript_id", "feature_id", "enstx", "isoform_id", "TXNAME"),
    required = TRUE,
    label = "transcript ID"
  )

  tx_id <- .fig6_strip_version(df[[transcript_col]])
  idx <- match(tx_id, annotation$transcripts$transcript_id)
  fields <- c(
    "transcript_id_full",
    "transcript_biotype",
    "gene_biotype",
    "gene_symbol",
    "transcript_name",
    "transcript_support_level"
  )
  for (nm in fields) {
    values <- annotation$transcripts[[nm]][idx]
    if (!nm %in% colnames(df)) {
      df[[nm]] <- values
    } else {
      replace <- is.na(df[[nm]]) | as.character(df[[nm]]) == ""
      df[[nm]][replace] <- values[replace]
    }
  }
  df$gencode_transcript_id <- tx_id
  df$gencode_gene_id <- annotation$transcripts$gene_id[idx]
  df
}

annotate_gene_biotypes <- function(df, annotation, gene_col = NULL) {
  if (is.null(df) || nrow(df) == 0) {
    return(df)
  }
  if (is.null(annotation$genes)) stop("Annotation object lacks $genes.")

  gene_col <- gene_col %||% .fig6_first_col(
    df,
    c("ensgene", "gene_id_full", "gene_id", "groupID", "GENEID"),
    required = TRUE,
    label = "gene ID"
  )
  gene_id <- .fig6_strip_version(df[[gene_col]])
  idx <- match(gene_id, annotation$genes$gene_id)

  if (!"gene_biotype" %in% colnames(df)) {
    df$gene_biotype <- annotation$genes$gene_biotype[idx]
  } else {
    replace <- is.na(df$gene_biotype) | as.character(df$gene_biotype) == ""
    df$gene_biotype[replace] <- annotation$genes$gene_biotype[idx][replace]
  }
  if (!"gene_symbol" %in% colnames(df)) {
    df$gene_symbol <- annotation$genes$gene_symbol[idx]
  }
  df$gencode_gene_id <- gene_id
  df
}

# --------------------------------------------------
# Locate and load cached analyses
# --------------------------------------------------

find_isoform_switch_cache <- function(
  switch_rds = NULL,
  prefer_strict = TRUE,
  required = TRUE
) {
  if (!is.null(switch_rds) && nzchar(switch_rds)) {
    return(.fig6_first_existing(switch_rds, required = required, label = "IsoformSwitchAnalyzeR cache"))
  }

  preferred_name <- if (isTRUE(prefer_strict)) {
    "switchList_filtered_strict_analyzed.rds"
  } else {
    "switchList_filtered_analyzed.rds"
  }
  alternate_name <- if (isTRUE(prefer_strict)) {
    "switchList_filtered_analyzed.rds"
  } else {
    "switchList_filtered_strict_analyzed.rds"
  }

  roots <- c(
    file.path(PROJECT_ROOT, "results", "analysis", "isoform_switch"),
    file.path(PROJECT_ROOT, "results", "analysis", "IsoformSwitchAnalyzeR"),
    file.path(PROJECT_ROOT, "results", "analysis", "isoformswitchanalyzer"),
    file.path(PROJECT_ROOT, "results", "IsoformSwitchAnalyzeR"),
    file.path(PROJECT_ROOT, "IsoformSwitchAnalyzeR"),
    file.path(PROJECT_ROOT, "resources", "cache")
  )
  candidates <- c(
    file.path(roots, preferred_name),
    file.path(roots, alternate_name)
  )

  if (!any(file.exists(candidates))) {
    found <- list.files(
      PROJECT_ROOT,
      pattern = "^switchList_filtered(_strict)?_analyzed\\.rds$",
      recursive = TRUE,
      full.names = TRUE
    )
    preferred <- found[basename(found) == preferred_name]
    candidates <- c(candidates, preferred, setdiff(found, preferred))
  }

  .fig6_first_existing(candidates, required = required, label = "IsoformSwitchAnalyzeR cache")
}

.fig6_find_deseq_cache <- function(deseq_rds = NULL, required = FALSE) {
  if (!is.null(deseq_rds) && nzchar(deseq_rds)) {
    return(.fig6_first_existing(deseq_rds, required = required, label = "DESeq2 cache"))
  }
  candidates <- c(
    file.path(PROJECT_ROOT, "results", "analysis", "deseq", "deseq_results.C1.rds"),
    Sys.glob(file.path(PROJECT_ROOT, "results", "analysis", "deseq", "deseq_results*.rds"))
  )
  .fig6_first_existing(candidates, required = required, label = "DESeq2 cache")
}

.fig6_standardize_features <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  tx_col <- .fig6_first_col(
    df,
    c("isoform_id", "transcript_id", "feature_id", "enstx"),
    required = TRUE,
    label = "IsoformSwitchAnalyzeR isoform ID"
  )
  gene_col <- .fig6_first_col(
    df,
    c("gene_id", "gene_name", "gene"),
    required = TRUE,
    label = "IsoformSwitchAnalyzeR gene ID"
  )
  c1_col <- .fig6_first_col(df, c("condition_1", "condition1"), required = TRUE)
  c2_col <- .fig6_first_col(df, c("condition_2", "condition2"), required = TRUE)

  if (tx_col != "isoform_id") df$isoform_id <- as.character(df[[tx_col]])
  if (gene_col != "gene_id") df$gene_id <- as.character(df[[gene_col]])
  if (c1_col != "condition_1") df$condition_1 <- as.character(df[[c1_col]])
  if (c2_col != "condition_2") df$condition_2 <- as.character(df[[c2_col]])

  df$isoform_id <- as.character(df$isoform_id)
  df$gene_id <- as.character(df$gene_id)
  df$condition_1 <- as.character(df$condition_1)
  df$condition_2 <- as.character(df$condition_2)
  df$comparison <- paste0(df$condition_2, "_vs_", df$condition_1)
  df$timepoint <- df$condition_2
  df
}

.fig6_subset_switch_object <- function(switch_list, keep) {
  .fig6_check_pkg("IsoformSwitchAnalyzeR")
  IsoformSwitchAnalyzeR::subsetSwitchAnalyzeRlist(switch_list, keep)
}

.fig6_extract_significant <- function(switch_list, condition_1, condition_2) {
  features <- .fig6_standardize_features(switch_list$isoformFeatures)
  keep <- features$condition_1 == condition_1 & features$condition_2 == condition_2
  if (!any(keep)) {
    return(data.frame())
  }

  sub_obj <- .fig6_subset_switch_object(switch_list, keep)
  top <- IsoformSwitchAnalyzeR::extractTopSwitches(
    sub_obj,
    extractGenes = FALSE,
    filterForConsequences = FALSE,
    n = NA,
    sortByQvals = TRUE
  )
  top <- as.data.frame(top, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(top) == 0) {
    return(data.frame())
  }

  top_tx_col <- .fig6_first_col(
    top,
    c("isoform_id", "transcript_id", "feature_id", "enstx"),
    required = TRUE,
    label = "extractTopSwitches isoform ID"
  )
  top$isoform_id <- as.character(top[[top_tx_col]])
  significant_ids <- unique(top$isoform_id)

  # isoformFeatures is the canonical source for comparison-specific dIF,
  # expression, and consequence columns. Supplement it with any columns returned
  # only by extractTopSwitches().
  out <- features[keep & features$isoform_id %in% significant_ids, , drop = FALSE]
  extra_cols <- setdiff(colnames(top), colnames(out))
  if (length(extra_cols) > 0) {
    top_extra <- top[, c("isoform_id", extra_cols), drop = FALSE]
    top_extra <- .fig6_make_unique(top_extra, "isoform_id")
    out <- merge(out, top_extra, by = "isoform_id", all.x = TRUE, sort = FALSE)
  }
  out$condition_1 <- condition_1
  out$condition_2 <- condition_2
  out$comparison <- paste0(condition_2, "_vs_", condition_1)
  out$timepoint <- condition_2
  out
}

.fig6_extract_all_significant <- function(switch_list) {
  features <- .fig6_standardize_features(switch_list$isoformFeatures)
  pairs <- unique(features[, c("condition_1", "condition_2"), drop = FALSE])
  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    .fig6_extract_significant(
      switch_list,
      condition_1 = pairs$condition_1[[i]],
      condition_2 = pairs$condition_2[[i]]
    )
  })
  .fig6_bind_rows(rows)
}

.fig6_native_isoform_annotations <- function(switch_list) {
  components <- names(switch_list)
  rows <- list()

  for (nm in components) {
    obj <- switch_list[[nm]]
    if (!is.data.frame(obj) || nrow(obj) == 0) next
    tx_col <- .fig6_first_col(
      obj,
      c("isoform_id", "transcript_id", "feature_id", "enstx"),
      required = FALSE
    )
    if (is.null(tx_col)) next

    useful <- grep(
      "nmd|ptc|premature|orf|coding|cpc|consequence",
      colnames(obj),
      ignore.case = TRUE,
      value = TRUE
    )
    useful <- setdiff(useful, tx_col)
    if (length(useful) == 0) next

    x <- as.data.frame(obj[, c(tx_col, useful), drop = FALSE], stringsAsFactors = FALSE)
    colnames(x)[1] <- "isoform_id"
    names(x)[-1] <- paste0(nm, "__", names(x)[-1])
    rows[[length(rows) + 1]] <- .fig6_make_unique(x, "isoform_id")
  }

  if (length(rows) == 0) {
    return(data.frame(isoform_id = character()))
  }
  out <- rows[[1]]
  if (length(rows) > 1) {
    for (i in 2:length(rows)) {
      out <- merge(out, rows[[i]], by = "isoform_id", all = TRUE, sort = FALSE)
    }
  }
  out
}

.fig6_attach_deseq <- function(switches, deseq_results, timepoints) {
  if (is.null(deseq_results) || nrow(switches) == 0) {
    switches$DEG <- NA
    return(switches)
  }

  full_by_tp <- deseq_results$results$gene_full_by_tp
  sig_by_tp <- deseq_results$results$gene_sig_by_tp
  rows <- lapply(timepoints, function(tp) {
    x <- switches[switches$condition_2 == tp, , drop = FALSE]
    if (nrow(x) == 0) {
      return(NULL)
    }

    full <- as.data.frame(full_by_tp[[tp]], stringsAsFactors = FALSE)
    sig <- as.data.frame(sig_by_tp[[tp]], stringsAsFactors = FALSE)
    id_values <- as.character(x$gene_id)

    full_gene_col <- .fig6_first_col(full, c("ensgene", "gene_id"), required = FALSE)
    full_symbol_col <- .fig6_first_col(full, c("symbol", "gene_name"), required = FALSE)
    lfc_col <- .fig6_first_col(
      full,
      c("log2FoldChange", "gene_log2_fold_change"),
      required = FALSE
    )

    idx <- rep(NA_integer_, nrow(x))
    if (!is.null(full_symbol_col)) {
      idx <- match(id_values, as.character(full[[full_symbol_col]]))
    }
    if (!is.null(full_gene_col)) {
      idx2 <- match(.fig6_strip_version(id_values), .fig6_strip_version(full[[full_gene_col]]))
      idx[is.na(idx)] <- idx2[is.na(idx)]
    }
    if (!is.null(lfc_col)) {
      existing <- if ("gene_log2_fold_change" %in% colnames(x)) {
        suppressWarnings(as.numeric(x$gene_log2_fold_change))
      } else {
        rep(NA_real_, nrow(x))
      }
      missing <- is.na(existing)
      existing[missing] <- suppressWarnings(as.numeric(full[[lfc_col]][idx]))[missing]
      x$gene_log2_fold_change <- existing
    }

    sig_keys <- character()
    sig_gene_col <- .fig6_first_col(sig, c("ensgene", "gene_id"), required = FALSE)
    sig_symbol_col <- .fig6_first_col(sig, c("symbol", "gene_name"), required = FALSE)
    if (!is.null(sig_gene_col)) {
      sig_keys <- c(sig_keys, .fig6_strip_version(sig[[sig_gene_col]]))
    }
    if (!is.null(sig_symbol_col)) sig_keys <- c(sig_keys, as.character(sig[[sig_symbol_col]]))
    x$DEG <- id_values %in% sig_keys | .fig6_strip_version(id_values) %in% sig_keys
    x
  })
  .fig6_bind_rows(rows)
}

.fig6_add_biotype_and_native_annotation <- function(switches, switch_list, gencode) {
  if (nrow(switches) == 0) {
    return(switches)
  }
  switches <- annotate_transcript_biotypes(
    switches,
    annotation = gencode,
    transcript_col = "isoform_id"
  )

  native <- .fig6_native_isoform_annotations(switch_list)
  if (nrow(native) > 0) {
    switches <- merge(switches, native, by = "isoform_id", all.x = TRUE, sort = FALSE)
  }

  if ("dIF" %in% colnames(switches)) {
    d_if <- suppressWarnings(as.numeric(switches$dIF))
    switches$usage_change_in_condition_2 <- ifelse(
      is.na(d_if),
      NA_character_,
      ifelse(d_if > 0, "increased", ifelse(d_if < 0, "decreased", "unchanged"))
    )
  } else {
    switches$usage_change_in_condition_2 <- NA_character_
  }
  switches$gencode_NMD_biotype <- switches$transcript_biotype == "nonsense_mediated_decay"
  switches
}

.fig6_biotype_enrichment <- function(
  significant,
  background,
  timepoints,
  reference_condition = "C1"
) {
  rows <- list()

  for (tp in timepoints) {
    sig_tp_all <- significant[
      significant$condition_1 == reference_condition & significant$condition_2 == tp, ,
      drop = FALSE
    ]
    bg_tp <- background[
      background$condition_1 == reference_condition & background$condition_2 == tp, ,
      drop = FALSE
    ]

    sig_tp_all <- .fig6_make_unique(sig_tp_all, "isoform_id")
    bg_tp <- .fig6_make_unique(bg_tp, "isoform_id")
    n_sig_total <- nrow(sig_tp_all)
    n_bg_total <- nrow(bg_tp)

    bg_mapped <- bg_tp[!is.na(bg_tp$transcript_biotype), , drop = FALSE]
    if (nrow(bg_mapped) == 0) next

    direction_sets <- list(
      all_switches = sig_tp_all,
      increased_in_hypoxia = sig_tp_all[
        !is.na(sig_tp_all$usage_change_in_condition_2) &
          sig_tp_all$usage_change_in_condition_2 == "increased", ,
        drop = FALSE
      ],
      decreased_in_hypoxia = sig_tp_all[
        !is.na(sig_tp_all$usage_change_in_condition_2) &
          sig_tp_all$usage_change_in_condition_2 == "decreased", ,
        drop = FALSE
      ]
    )

    one_direction <- lapply(names(direction_sets), function(direction_name) {
      sig_unmapped_included <- direction_sets[[direction_name]]
      sig_mapped <- sig_unmapped_included[
        !is.na(sig_unmapped_included$transcript_biotype), ,
        drop = FALSE
      ]
      sig_total <- nrow(sig_mapped)
      bg_total <- nrow(bg_mapped)
      biotypes <- sort(unique(bg_mapped$transcript_biotype))

      one_test <- lapply(biotypes, function(bt) {
        sig_in <- sum(sig_mapped$transcript_biotype == bt)
        bg_in <- sum(bg_mapped$transcript_biotype == bt)

        a <- sig_in
        b <- sig_total - sig_in
        c <- max(bg_in - sig_in, 0)
        d <- max((bg_total - bg_in) - b, 0)
        ft <- if (sig_total > 0) {
          stats::fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
        } else {
          NULL
        }

        data.frame(
          timepoint = tp,
          comparison = paste0(tp, "_vs_", reference_condition),
          switch_direction = direction_name,
          transcript_biotype = bt,
          n_switching_transcripts = sig_in,
          n_switching_transcripts_with_biotype = sig_total,
          n_total_switching_transcripts = nrow(sig_unmapped_included),
          switching_biotype_mapping_fraction = if (nrow(sig_unmapped_included) > 0) {
            sig_total / nrow(sig_unmapped_included)
          } else {
            NA_real_
          },
          fraction_switching = if (sig_total > 0) sig_in / sig_total else NA_real_,
          n_background_transcripts = bg_in,
          n_background_transcripts_with_biotype = bg_total,
          n_total_background_transcripts = n_bg_total,
          background_biotype_mapping_fraction = bg_total / n_bg_total,
          fraction_background = bg_in / bg_total,
          odds_ratio = if (is.null(ft)) NA_real_ else unname(ft$estimate),
          pvalue = if (is.null(ft)) NA_real_ else ft$p.value,
          stringsAsFactors = FALSE
        )
      })
      one_test <- do.call(rbind, one_test)
      one_test$padj <- stats::p.adjust(one_test$pvalue, method = "BH")
      one_test$enriched_at_fdr_0.10 <- !is.na(one_test$padj) &
        one_test$padj < 0.10 &
        one_test$odds_ratio > 1
      one_test$depleted_at_fdr_0.10 <- !is.na(one_test$padj) &
        one_test$padj < 0.10 &
        one_test$odds_ratio < 1
      one_test
    })

    rows[[tp]] <- .fig6_bind_rows(one_direction)
  }

  .fig6_bind_rows(rows)
}

.fig6_consequence_summary <- function(
  switch_list,
  alpha = 0.10,
  reference_condition = "C1",
  timepoints = c("H1", "H3", "H24")
) {
  consequences <- c(
    "5_utr_seq_similarity",
    "3_utr_seq_similarity",
    "isoform_seq_similarity",
    "domains_identified",
    "domain_isotype",
    "exon_number",
    "IDR_identified",
    "intron_retention",
    "isoform_topology",
    "extracellular_region_count",
    "intracellular_region_count",
    "last_exon",
    "sub_cell_location",
    "sub_cell_shift_to_cytoplasm",
    "sub_cell_shift_to_nucleus",
    "tss",
    "tts"
  )

  features <- .fig6_standardize_features(switch_list$isoformFeatures)
  keep <- features$condition_1 == reference_condition &
    features$condition_2 %in% timepoints
  sub_obj <- .fig6_subset_switch_object(switch_list, keep)

  summary_plot <- IsoformSwitchAnalyzeR::extractConsequenceSummary(
    sub_obj,
    consequencesToAnalyze = consequences,
    includeCombined = FALSE,
    plotGenes = TRUE,
    removeEmptyConsequences = TRUE,
    asFractionTotal = TRUE,
    returnResult = FALSE
  )
  enrichment <- IsoformSwitchAnalyzeR::extractConsequenceEnrichment(
    sub_obj,
    consequencesToAnalyze = consequences,
    alpha = alpha,
    analysisOppositeConsequence = TRUE,
    returnResult = TRUE
  )

  dat <- as.data.frame(summary_plot$data, stringsAsFactors = FALSE, check.names = FALSE)
  dat$plotComparison <- gsub("[\r\n]+", " ", dat$plotComparison)
  for (tp in timepoints) {
    dat$plotComparison <- gsub(
      paste0(reference_condition, " vs ", tp),
      tp,
      dat$plotComparison,
      fixed = TRUE
    )
  }
  dat <- dat[
    !grepl("Mixed domain isotype changes|Domain switch|IDR switch", dat$switchConsequence), ,
    drop = FALSE
  ]

  enrichment <- as.data.frame(enrichment, stringsAsFactors = FALSE, check.names = FALSE)
  enrichment$feature <- gsub("\\s*\\([^\\)]+\\)", "", enrichment$feature)
  keep_cols <- intersect(c("condition_2", "feature", "Significant"), colnames(enrichment))
  enrichment <- enrichment[, keep_cols, drop = FALSE]
  if (all(c("condition_2", "feature", "Significant") %in% colnames(enrichment))) {
    enrichment <- enrichment[!duplicated(enrichment[, c("condition_2", "feature")]), , drop = FALSE]
    dat <- merge(
      dat,
      enrichment,
      by.x = c("plotComparison", "switchConsequence"),
      by.y = c("condition_2", "feature"),
      all.x = TRUE,
      sort = FALSE
    )
  } else {
    dat$Significant <- FALSE
  }
  dat$Significant[is.na(dat$Significant)] <- FALSE
  dat$Significant <- as.logical(dat$Significant)

  list(
    data = dat,
    mapping = summary_plot$mapping,
    enrichment = enrichment,
    consequences = consequences
  )
}

run_isoform_switch_analysis <- function(
  switch_list = NULL,
  switch_rds = NULL,
  deseq_results = NULL,
  deseq_rds = NULL,
  gtf_path = NULL,
  reference_condition = "C1",
  timepoints = c("H1", "H3", "H24"),
  consequence_alpha = 0.10,
  outdir = NULL,
  force = FALSE,
  save_tables = TRUE
) {
  .fig6_check_pkg("IsoformSwitchAnalyzeR")
  outdir <- .fig6_outdir(outdir)
  analysis_rds <- file.path(outdir, "fig_6_isoform_switch_analysis.rds")

  switch_rds <- if (is.null(switch_list)) {
    find_isoform_switch_cache(switch_rds = switch_rds, prefer_strict = TRUE)
  } else {
    switch_rds %||% NA_character_
  }
  if (is.null(switch_list)) {
    message("[Figure 6] Loading IsoformSwitchAnalyzeR cache: ", switch_rds)
    switch_list <- readRDS(switch_rds)
  }
  if (is.null(switch_list$isoformFeatures)) {
    stop("The switch object does not contain $isoformFeatures.")
  }

  if (is.null(deseq_results)) {
    deseq_path <- .fig6_find_deseq_cache(deseq_rds, required = FALSE)
    if (!is.na(deseq_path)) {
      message("[Figure 6] Loading DESeq2 cache: ", deseq_path)
      deseq_results <- readRDS(deseq_path)
    }
  }

  can_use_cache <- !is.na(switch_rds) && nzchar(switch_rds)
  if (!isTRUE(force) && can_use_cache && file.exists(analysis_rds)) {
    cached <- readRDS(analysis_rds)
    if (
      identical(cached$meta$analysis_version, .fig6_analysis_version) &&
        identical(cached$meta$timepoints, timepoints) &&
        identical(cached$meta$reference_condition, reference_condition) &&
        identical(
          normalizePath(cached$meta$switch_rds, mustWork = FALSE),
          normalizePath(switch_rds, mustWork = FALSE)
        )
    ) {
      message("[Figure 6] Loading cached analysis: ", analysis_rds)
      cached$switch_list <- switch_list
      cached$deseq_results <- deseq_results
      return(cached)
    }
  }

  gencode <- get_gencode_biotype_annotation(gtf_path = gtf_path)
  background <- .fig6_standardize_features(switch_list$isoformFeatures)
  background <- background[
    background$condition_1 == reference_condition &
      background$condition_2 %in% timepoints, ,
    drop = FALSE
  ]
  background <- annotate_transcript_biotypes(
    background,
    annotation = gencode,
    transcript_col = "isoform_id"
  )

  message("[Figure 6] Extracting significant isoform switches")
  all_switches <- .fig6_extract_all_significant(switch_list)
  c1_switches <- all_switches[
    all_switches$condition_1 == reference_condition &
      all_switches$condition_2 %in% timepoints, ,
    drop = FALSE
  ]
  c1_switches <- .fig6_add_biotype_and_native_annotation(
    c1_switches,
    switch_list = switch_list,
    gencode = gencode
  )
  c1_switches <- .fig6_attach_deseq(c1_switches, deseq_results, timepoints)

  biotype_enrichment <- .fig6_biotype_enrichment(
    significant = c1_switches,
    background = background,
    timepoints = timepoints,
    reference_condition = reference_condition
  )

  consequence <- tryCatch(
    .fig6_consequence_summary(
      switch_list,
      alpha = consequence_alpha,
      reference_condition = reference_condition,
      timepoints = timepoints
    ),
    error = function(e) {
      warning(
        "Could not recreate the consequence summary from the cached object: ",
        conditionMessage(e)
      )
      list(
        data = data.frame(),
        mapping = NULL,
        enrichment = data.frame(),
        consequences = character(),
        error = conditionMessage(e)
      )
    }
  )

  summary_tbl <- do.call(rbind, lapply(timepoints, function(tp) {
    x <- c1_switches[c1_switches$condition_2 == tp, , drop = FALSE]
    x <- .fig6_make_unique(x, "isoform_id")
    mapped <- !is.na(x$transcript_biotype) & x$transcript_biotype != ""
    n_mapped <- sum(mapped)
    n_nmd <- sum(x$gencode_NMD_biotype, na.rm = TRUE)
    data.frame(
      timepoint = tp,
      comparison = paste0(tp, "_vs_", reference_condition),
      significant_switching_transcripts = length(unique(x$isoform_id)),
      significant_switching_genes = length(unique(x$gene_id)),
      transcripts_with_gencode_biotype = n_mapped,
      transcript_biotype_mapping_fraction = if (nrow(x) > 0) n_mapped / nrow(x) else NA_real_,
      nmd_biotype_transcripts = n_nmd,
      nmd_biotype_fraction_among_mapped = if (n_mapped > 0) n_nmd / n_mapped else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  out <- list(
    meta = list(
      analysis_version = .fig6_analysis_version,
      switch_rds = switch_rds,
      reference_condition = reference_condition,
      timepoints = timepoints,
      consequence_alpha = consequence_alpha,
      gencode_gtf = gencode$meta$source$path,
      biotype_background = paste(
        "All unique transcripts in isoformFeatures eligible for each",
        paste0(reference_condition, "-to-hypoxia contrast")
      )
    ),
    tables = list(
      isoform_switches = c1_switches,
      biotype_enrichment = biotype_enrichment,
      background_transcripts = background,
      all_pairwise_switches = all_switches,
      consequence_summary = consequence$data,
      consequence_enrichment = consequence$enrichment
    ),
    summary = summary_tbl,
    consequence = consequence,
    deseq_results = deseq_results,
    switch_list = switch_list
  )

  cache_out <- out
  cache_out$switch_list <- NULL
  cache_out$deseq_results <- NULL
  saveRDS(cache_out, analysis_rds)

  if (isTRUE(save_tables)) {
    utils::write.table(
      c1_switches,
      file = file.path(outdir, "fig_6_isoform_switches.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    utils::write.table(
      biotype_enrichment,
      file = file.path(outdir, "fig_6_transcript_biotype_enrichment.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    utils::write.table(
      summary_tbl,
      file = file.path(outdir, "fig_6_isoform_switch_summary.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  out
}

# --------------------------------------------------
# Figure 6A - Switch overlap and DEG fraction
# --------------------------------------------------

.fig6_pair_gene_sets <- function(all_switches) {
  pairs <- list(
    "C1 vs H1" = c("C1", "H1"),
    "C1 vs H3" = c("C1", "H3"),
    "C1 vs H24" = c("C1", "H24"),
    "H1 vs H3" = c("H1", "H3"),
    "H1 vs H24" = c("H1", "H24"),
    "H3 vs H24" = c("H3", "H24")
  )
  lapply(pairs, function(pair) {
    keep <- (
      all_switches$condition_1 == pair[[1]] &
        all_switches$condition_2 == pair[[2]]
    ) | (
      all_switches$condition_1 == pair[[2]] &
        all_switches$condition_2 == pair[[1]]
    )
    unique(as.character(all_switches$gene_id[keep]))
  })
}

.fig6_deg_keys <- function(deseq_results, timepoints) {
  if (is.null(deseq_results)) {
    return(character())
  }
  sig_by_tp <- deseq_results$results$gene_sig_by_tp
  unique(unlist(lapply(timepoints, function(tp) {
    x <- as.data.frame(sig_by_tp[[tp]], stringsAsFactors = FALSE)
    gene_col <- .fig6_first_col(x, c("ensgene", "gene_id"), required = FALSE)
    symbol_col <- .fig6_first_col(x, c("symbol", "gene_name"), required = FALSE)
    c(
      if (is.null(gene_col)) character() else .fig6_strip_version(x[[gene_col]]),
      if (is.null(symbol_col)) character() else as.character(x[[symbol_col]])
    )
  }), use.names = FALSE))
}

plot_switch_overlap_and_deg_fraction <- function(
  res,
  base_size = 12,
  outdir = NULL
) {
  .fig6_check_pkg(c("ggplot2", "ComplexUpset", "patchwork", "scales"))
  outdir <- .fig6_outdir(outdir)

  sets <- .fig6_pair_gene_sets(res$tables$all_pairwise_switches)
  all_genes <- sort(unique(unlist(sets, use.names = FALSE)))
  upset_df <- data.frame(gene_id = all_genes, stringsAsFactors = FALSE)
  for (nm in names(sets)) upset_df[[nm]] <- all_genes %in% sets[[nm]]

  set_colors <- c(
    "C1 vs H1" = "#377EB8",
    "C1 vs H3" = "#66A61E",
    "C1 vs H24" = "#7570B3",
    "H1 vs H3" = "#A6A6A6",
    "H1 vs H24" = "#737373",
    "H3 vs H24" = "#404040"
  )

  p_upset <- ComplexUpset::upset(
    upset_df,
    intersect = names(sets),
    name = "Contrasts",
    width_ratio = 0.22,
    set_sizes = ComplexUpset::upset_set_size() +
      ggplot2::ylab("Set size"),
    queries = lapply(names(sets), function(nm) {
      ComplexUpset::upset_query(set = nm, fill = unname(set_colors[[nm]]))
    }),
    themes = ComplexUpset::upset_default_themes(
      text = ggplot2::element_text(size = base_size)
    )
  )

  c1_genes <- unique(res$tables$isoform_switches$gene_id)
  deg_keys <- .fig6_deg_keys(res$deseq_results, res$meta$timepoints)
  is_deg <- c1_genes %in% deg_keys | .fig6_strip_version(c1_genes) %in% deg_keys
  pie <- data.frame(
    group = c("Switch only", "Differentially expressed"),
    n = c(sum(!is_deg), sum(is_deg)),
    stringsAsFactors = FALSE
  )
  pie$proportion <- pie$n / sum(pie$n)
  pie$label <- paste0(pie$group, "\n", scales::percent(pie$proportion, accuracy = 1))

  p_pie <- ggplot2::ggplot(pie, ggplot2::aes(x = "", y = n, fill = group)) +
    ggplot2::geom_col(color = "white", linewidth = 0.4) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      position = ggplot2::position_stack(vjust = 0.5),
      color = "white",
      fontface = "bold",
      size = 3.4
    ) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::scale_fill_manual(
      values = c("Switch only" = "#C84655", "Differentially expressed" = "#202020")
    ) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(legend.position = "none")

  p <- patchwork::wrap_plots(p_upset, p_pie, widths = c(3.1, 1))
  ggplot2::ggsave(
    file.path(outdir, "fig_6_A_switch_overlap_and_deg_fraction.pdf"),
    p,
    width = 10.5,
    height = 5.0
  )
  p
}

# --------------------------------------------------
# Figure 6B - Example switches
# --------------------------------------------------

.fig6_placeholder <- function(label, base_size = 12) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = label, size = 4) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::theme_void(base_size = base_size)
}

plot_isoform_switch_examples <- function(
  res,
  example_gene = "PFKFB3",
  example_isoforms = c("ENST00000676698.1", "ENST00000519627.2"),
  condition_1 = "C1",
  condition_2 = "H3",
  base_size = 10,
  outdir = NULL
) {
  .fig6_check_pkg(c("IsoformSwitchAnalyzeR", "ggplot2", "patchwork"))
  outdir <- .fig6_outdir(outdir)
  switch_list <- res$switch_list
  if (is.null(switch_list)) stop("The analysis object does not include the loaded switch list.")

  p_gene <- tryCatch(
    IsoformSwitchAnalyzeR::switchPlot(
      switch_list,
      gene = example_gene,
      condition1 = condition_1,
      condition2 = condition_2,
      localTheme = ggplot2::theme_bw(base_size = base_size)
    ),
    error = function(e) {
      warning("Could not plot ", example_gene, ": ", conditionMessage(e))
      .fig6_placeholder(paste(example_gene, conditionMessage(e), sep = "\n"), base_size)
    }
  )

  p_isoforms <- tryCatch(
    IsoformSwitchAnalyzeR::switchPlot(
      switch_list,
      isoform_id = example_isoforms,
      condition1 = condition_1,
      condition2 = condition_2,
      localTheme = ggplot2::theme_bw(base_size = base_size)
    ),
    error = function(e) {
      warning("Could not plot requested EIF3E isoforms: ", conditionMessage(e))
      .fig6_placeholder(paste("EIF3E", conditionMessage(e), sep = "\n"), base_size)
    }
  )

  p <- patchwork::wrap_plots(p_gene, p_isoforms, ncol = 1)
  ggplot2::ggsave(
    file.path(outdir, "fig_6_B_isoform_switch_examples.pdf"),
    p,
    width = 8.0,
    height = 8.4
  )
  p
}

# --------------------------------------------------
# Figure 6C - Switch usage versus gene expression
# --------------------------------------------------

plot_switch_usage_vs_gene_expression <- function(
  res,
  base_size = 12,
  outdir = NULL
) {
  .fig6_check_pkg("ggplot2")
  outdir <- .fig6_outdir(outdir)
  df <- as.data.frame(res$tables$isoform_switches, stringsAsFactors = FALSE)

  if (!all(c("gene_log2_fold_change", "dIF") %in% colnames(df))) {
    stop("Figure 6C requires gene_log2_fold_change and dIF columns.")
  }
  df$gene_log2_fold_change <- suppressWarnings(as.numeric(df$gene_log2_fold_change))
  df$dIF <- suppressWarnings(as.numeric(df$dIF))
  df$DEG <- factor(as.logical(df$DEG), levels = c(FALSE, TRUE))
  df$condition_2 <- factor(df$condition_2, levels = res$meta$timepoints)

  tp_labels <- c(H1 = "1 hour", H3 = "3 hours", H24 = "24 hours")
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = gene_log2_fold_change, y = dIF, color = DEG)
  ) +
    ggplot2::geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey55") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    ggplot2::geom_point(alpha = 0.75, size = 1.3) +
    ggplot2::facet_wrap(
      ~condition_2,
      nrow = 1,
      labeller = ggplot2::as_labeller(tp_labels)
    ) +
    ggplot2::scale_color_manual(
      name = "Differentially\nexpressed gene",
      values = c("FALSE" = "#C84655", "TRUE" = "#202020"),
      labels = c("FALSE" = "No", "TRUE" = "Yes"),
      na.value = "grey70",
      drop = FALSE
    ) +
    ggplot2::labs(
      x = expression("Gene expression (log"[2] * " fold change)"),
      y = "Difference in isoform usage (dIF)"
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "white"),
      legend.position = "right"
    )

  ggplot2::ggsave(
    file.path(outdir, "fig_6_C_switch_usage_vs_gene_expression.pdf"),
    p,
    width = 8.0,
    height = 3.4
  )
  p
}

# --------------------------------------------------
# Figure 6D - Consequence summary
# --------------------------------------------------

plot_switch_consequence_summary <- function(
  res,
  base_size = 9,
  outdir = NULL
) {
  .fig6_check_pkg(c("ggplot2", "scales"))
  outdir <- .fig6_outdir(outdir)
  dat <- as.data.frame(res$consequence$data, stringsAsFactors = FALSE, check.names = FALSE)
  mapping <- res$consequence$mapping

  if (nrow(dat) == 0 || is.null(mapping)) {
    p <- .fig6_placeholder(
      "Consequence summary unavailable in the cached switch object.",
      base_size = base_size
    )
  } else {
    feature_order <- c(
      "Tss",
      "Exon\nnumber",
      "Intron\nretention",
      "Last\nexon",
      "Tts",
      "Isoform seq\nsimilarity",
      "5 utr seq\nsimilarity",
      "3 utr seq\nsimilarity",
      "Domains\nidentified",
      "Domain\nisotype",
      "IDR\nidentified",
      "Isoform\ntopology",
      "Intracellular\nregion count",
      "Extracellular\nregion count",
      "Sub cell\nlocation",
      "Sub cell shift\nto cytoplasm",
      "Sub cell shift\nto nucleus"
    )

    p <- ggplot2::ggplot(dat, mapping) +
      ggplot2::geom_col(ggplot2::aes(fill = Significant)) +
      ggplot2::scale_fill_manual(
        "Significant enrichment",
        values = c("TRUE" = "#C84655", "FALSE" = "#202020"),
        breaks = c(TRUE, FALSE),
        labels = c("Yes", "No")
      ) +
      ggplot2::scale_y_continuous(labels = scales::percent) +
      ggplot2::facet_grid(
        factor(plotComparison, levels = res$meta$timepoints) ~
          factor(featureCompared, levels = feature_order),
        scales = "free",
        space = "free_x"
      ) +
      ggplot2::labs(
        x = "Consequence of switch (feature of the hypoxia-favored isoform)",
        y = "Fraction of switched genes"
      ) +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        strip.text.y = ggplot2::element_text(angle = 0),
        strip.background = ggplot2::element_rect(fill = "white"),
        axis.text.x = ggplot2::element_text(angle = -45, hjust = 0, vjust = 1),
        legend.position = "bottom"
      )
  }

  ggplot2::ggsave(
    file.path(outdir, "fig_6_D_switch_consequence_summary.pdf"),
    p,
    width = 11.0,
    height = 5.6
  )
  p
}

# --------------------------------------------------
# Convenience wrapper
# --------------------------------------------------

plot_fig_6_all <- function(
  switch_list = NULL,
  switch_rds = NULL,
  deseq_results = NULL,
  deseq_rds = NULL,
  gtf_path = NULL,
  reference_condition = "C1",
  timepoints = c("H1", "H3", "H24"),
  consequence_alpha = 0.10,
  example_gene = "PFKFB3",
  example_isoforms = c("ENST00000676698.1", "ENST00000519627.2"),
  outdir = NULL,
  force = FALSE,
  save_tables = TRUE
) {
  .fig6_check_pkg(c("ggplot2", "patchwork"))
  outdir <- .fig6_outdir(outdir)

  res <- run_isoform_switch_analysis(
    switch_list = switch_list,
    switch_rds = switch_rds,
    deseq_results = deseq_results,
    deseq_rds = deseq_rds,
    gtf_path = gtf_path,
    reference_condition = reference_condition,
    timepoints = timepoints,
    consequence_alpha = consequence_alpha,
    outdir = outdir,
    force = force,
    save_tables = save_tables
  )

  p6a <- plot_switch_overlap_and_deg_fraction(res, outdir = outdir)
  p6b <- plot_isoform_switch_examples(
    res,
    example_gene = example_gene,
    example_isoforms = example_isoforms,
    condition_1 = reference_condition,
    condition_2 = "H3",
    outdir = outdir
  )
  p6c <- plot_switch_usage_vs_gene_expression(res, outdir = outdir)
  p6d <- plot_switch_consequence_summary(res, outdir = outdir)

  design <- "
  ABB
  CBB
  DDD
  DDD
  "
  combined <- patchwork::wrap_plots(
    A = patchwork::wrap_elements(full = p6a),
    B = patchwork::wrap_elements(full = p6b),
    C = p6c,
    D = p6d,
    design = design
  ) +
    patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 14))

  ggplot2::ggsave(
    file.path(outdir, "fig_6.pdf"),
    combined,
    width = 8.5,
    height = 11.25
  )

  invisible(list(
    analysis = res,
    fig_6_a = p6a,
    fig_6_b = p6b,
    fig_6_c = p6c,
    fig_6_d = p6d,
    combined = combined
  ))
}
