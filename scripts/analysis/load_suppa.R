PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
if (!nzchar(PROJECT_ROOT)) stop("PROJECT_ROOT env var not set.")

HELPERS_FILE <- file.path(PROJECT_ROOT, "scripts/utils/helpers.R")
if (!file.exists(HELPERS_FILE)) stop("helpers.R not found at: ", HELPERS_FILE)
source(HELPERS_FILE)

# ------------------------------------------------------------------------------
# Internal Helpers
# ------------------------------------------------------------------------------

.parse_suppa_event_id <- function(event_id) {
  event_id <- as.character(event_id)

  ensgene_raw <- sub(";.*$", "", event_id)
  event_part <- sub("^[^;]+;", "", event_id)

  parts <- strsplit(event_part, ":", fixed = TRUE)
  n <- length(parts)

  event_type <- rep(NA_character_, n)
  chr <- rep(NA_character_, n)
  coord1 <- rep(NA_character_, n)
  coord2 <- rep(NA_character_, n)
  coord3 <- rep(NA_character_, n)
  coord4 <- rep(NA_character_, n)
  strand <- rep(NA_character_, n)

  for (i in seq_len(n)) {
    x <- parts[[i]]
    nx <- length(x)

    if (nx >= 1) event_type[i] <- x[1]
    if (nx >= 2) chr[i] <- x[2]
    if (nx >= 3) strand[i] <- x[nx]

    if (nx >= 4) {
      coords <- x[3:(nx - 1)]
      if (length(coords) >= 1) coord1[i] <- coords[1]
      if (length(coords) >= 2) coord2[i] <- coords[2]
      if (length(coords) >= 3) coord3[i] <- coords[3]
      if (length(coords) >= 4) coord4[i] <- coords[4]
    }
  }

  data.frame(
    ensgene_raw = ensgene_raw,
    event_type = event_type,
    chr = chr,
    coord1 = coord1,
    coord2 = coord2,
    coord3 = coord3,
    coord4 = coord4,
    strand = strand,
    stringsAsFactors = FALSE
  )
}

.read_one_suppa_dpsi <- function(path) {
  if (!file.exists(path)) stop("SUPPA dpsi file not found: ", path)

  # SUPPA dpsi files often have:
  # header = dPSI / p-value columns only
  # rows   = event_id + values
  x <- read.delim(
    path,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (ncol(x) < 2) {
    stop("Expected at least 2 columns in SUPPA dpsi file: ", path)
  }

  out <- as.data.frame(x, stringsAsFactors = FALSE)
  out$event_id <- rownames(out)
  rownames(out) <- NULL

  cn <- colnames(out)

  dpsi_col <- grep("dpsi", cn, ignore.case = TRUE, value = TRUE)
  pval_col <- grep("p[-._]?val|pvalue", cn, ignore.case = TRUE, value = TRUE)

  if (length(dpsi_col) == 0) dpsi_col <- cn[1]
  if (length(pval_col) == 0) pval_col <- cn[2]

  dpsi_col <- dpsi_col[1]
  pval_col <- pval_col[1]

  colnames(out)[match(dpsi_col, colnames(out))] <- "dpsi"
  colnames(out)[match(pval_col, colnames(out))] <- "pvalue"

  out$dpsi <- suppressWarnings(as.numeric(out$dpsi))
  out$pvalue <- suppressWarnings(as.numeric(out$pvalue))

  parsed <- .parse_suppa_event_id(out$event_id)

  out <- data.frame(
    event_id = out$event_id,
    parsed,
    dpsi = out$dpsi,
    pvalue = out$pvalue,
    stringsAsFactors = FALSE
  )

  out$ensgene <- .strip_ens_version(out$ensgene_raw)

  out
}

.read_suppa_isotpm <- function(path) {
  x <- utils::read.delim(
    path,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  sample_names <- colnames(x)

  # Because SUPPA files have sample names only in the header,
  # read.delim treats the first transcript column as rownames sometimes.
  if (!grepl("^ENS", x[[1]][1])) {
    x <- tibble::rownames_to_column(x, var = "transcript_id")
  }

  list(
    df = x,
    sample_names = sample_names
  )
}

.write_suppa_isotpm <- function(df, path, sample_names = NULL) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  df <- as.data.frame(df, check.names = FALSE)

  if (!is.null(sample_names)) {
    stopifnot(length(sample_names) == ncol(df) - 1)
  } else {
    sample_names <- colnames(df)[-1]
  }

  # First column should be transcript IDs
  tx <- df[[1]]
  mat <- df[, -1, drop = FALSE]

  # Write SUPPA-style file:
  # header = sample names only
  # rows   = transcript ID + values
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)

  writeLines(paste(sample_names, collapse = "\t"), con = con)

  utils::write.table(
    data.frame(tx, mat, check.names = FALSE),
    file = con,
    quote = FALSE,
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE
  )

  invisible(path)
}

.filter_suppa_isotpm_by_tx <- function(isotpm_df, tx_universe) {
  tx_universe <- unique(.strip_version(tx_universe))

  tx_ids <- .strip_version(isotpm_df[[1]])
  isotpm_df[tx_ids %in% tx_universe, , drop = FALSE]
}

.standardize_suppa_res <- function(res_df,
                                   timepoint,
                                   contrast,
                                   tool = "suppa2",
                                   annot_df = NULL,
                                   p_cutoff = 0.10) {
  out <- res_df

  out$timepoint <- timepoint
  out$contrast <- contrast
  out$tool <- tool

  out <- .attach_gene_annot(out, annot_df = annot_df)

  out$direction <- ifelse(
    is.na(out$pvalue), NA_character_,
    ifelse(out$pvalue < p_cutoff & out$dpsi > 0, "up",
      ifelse(out$pvalue < p_cutoff & out$dpsi < 0, "down", "ns")
    )
  )

  out$neglog10_pvalue <- ifelse(is.na(out$pvalue), NA_real_, -log10(out$pvalue))

  # convenience aliases to parallel DESeq2 output style
  out$gene_id <- out$ensgene
  out$comparison <- contrast

  out
}

.run_one_suppa_timepoint <- function(suppa_dir,
                                     timepoint,
                                     ref_level,
                                     p_cutoff = 0.10,
                                     annot_df = NULL) {
  dpsi_path <- file.path(suppa_dir, "diffSplice", paste0("diffSplice_", timepoint, ".dpsi"))
  contrast <- paste0(timepoint, "_vs_", ref_level)

  raw <- .read_one_suppa_dpsi(dpsi_path)

  full <- .standardize_suppa_res(
    res_df = raw,
    timepoint = timepoint,
    contrast = contrast,
    annot_df = annot_df,
    p_cutoff = p_cutoff
  )

  sig <- full[!is.na(full$pvalue) & full$pvalue < p_cutoff, , drop = FALSE]
  sig <- sig[order(sig$pvalue, -abs(sig$dpsi)), , drop = FALSE]

  sig_genes <- unique(sig[, c("ensgene", "timepoint", "contrast", "tool"), drop = FALSE])
  sig_genes <- sig_genes[order(sig_genes$ensgene), , drop = FALSE]

  list(
    contrast = contrast,
    full = full,
    sig = sig,
    sig_genes = sig_genes
  )
}

.make_suppa_summary <- function(sig_event_list,
                                sig_gene_list,
                                p_cutoff = 0.10) {
  tps <- names(sig_event_list)

  do.call(rbind, lapply(tps, function(tp) {
    ev <- sig_event_list[[tp]]
    gn <- sig_gene_list[[tp]]

    data.frame(
      timepoint = tp,
      p_cutoff = p_cutoff,
      n_sig_events = nrow(ev),
      n_sig_genes = nrow(gn),
      n_up_events = sum(ev$direction == "up", na.rm = TRUE),
      n_down_events = sum(ev$direction == "down", na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

.write_suppa_xlsx <- function(sig_event_list,
                              sig_gene_list,
                              summary_tbl,
                              out_xlsx) {
  .check_pkg("openxlsx")
  wb <- openxlsx::createWorkbook()

  openxlsx::addWorksheet(wb, "Summary")
  openxlsx::writeData(wb, "Summary", summary_tbl)

  for (tp in names(sig_event_list)) {
    ev_sheet <- paste0(tp, "_events")
    gn_sheet <- paste0(tp, "_genes")

    openxlsx::addWorksheet(wb, ev_sheet)
    if (nrow(sig_event_list[[tp]]) == 0) {
      openxlsx::writeData(wb, ev_sheet, data.frame(note = "No significant events at this cutoff."))
    } else {
      openxlsx::writeData(wb, ev_sheet, sig_event_list[[tp]])
    }

    openxlsx::addWorksheet(wb, gn_sheet)
    if (nrow(sig_gene_list[[tp]]) == 0) {
      openxlsx::writeData(wb, gn_sheet, data.frame(note = "No significant genes at this cutoff."))
    } else {
      openxlsx::writeData(wb, gn_sheet, sig_gene_list[[tp]])
    }
  }

  dir.create(dirname(out_xlsx), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
  invisible(TRUE)
}

.summarize_splicing_by_event_type <- function(suppa_res,
                                              sig_only = TRUE) {
  if (!("event_type" %in% colnames(suppa_res$results$full_all))) {
    stop("event_type column not found in SUPPA results.")
  }

  df <- if (sig_only) {
    suppa_res$results$sig_all
  } else {
    suppa_res$results$full_all
  }

  out <- aggregate(
    list(n_events = df$event_type),
    by = list(
      timepoint = df$timepoint,
      event_type = df$event_type
    ),
    FUN = length
  )

  out <- out[order(out$timepoint, out$event_type), ]
  rownames(out) <- NULL

  out
}

# ------------------------------------------------------------------------------
# Public API
# ------------------------------------------------------------------------------

load_suppa <- function(timepoints = c("H1", "H3", "H24"),
                       ref_level = "C1",
                       p_cutoff = 0.10,
                       annot_df = NULL,
                       force = FALSE,
                       suppa_dir = NULL,
                       out_rds = NULL,
                       out_xlsx = NULL) {
  if (is.null(suppa_dir)) {
    suppa_dir <- file.path(
      .get_results_dir(),
      "analysis/suppa/standard",
      ref_level
    )
  }

  if (!dir.exists(suppa_dir)) {
    stop("SUPPA directory not found: ", suppa_dir)
  }

  outdir <- suppa_dir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(out_rds)) out_rds <- file.path(outdir, "suppa_results.rds")
  if (is.null(out_xlsx)) out_xlsx <- file.path(outdir, "suppa_results.xlsx")

  if (file.exists(out_rds) && !force) {
    message("Loading cached SUPPA2 results: ", out_rds)
    return(readRDS(out_rds))
  }

  if (is.null(annot_df)) {
    annot_df <- .get_annot()
  }

  by_tp <- setNames(lapply(timepoints, function(tp) {
    .run_one_suppa_timepoint(
      suppa_dir = suppa_dir,
      timepoint = tp,
      ref_level = ref_level,
      p_cutoff = p_cutoff,
      annot_df = annot_df
    )
  }), timepoints)

  full_list <- lapply(by_tp, `[[`, "full")
  sig_event_list <- lapply(by_tp, `[[`, "sig")
  sig_gene_list <- lapply(by_tp, `[[`, "sig_genes")

  full_all <- do.call(rbind, full_list)
  sig_all <- do.call(rbind, sig_event_list)
  sig_genes_all <- unique(do.call(rbind, sig_gene_list))

  summary_tbl <- .make_suppa_summary(
    sig_event_list = sig_event_list,
    sig_gene_list = sig_gene_list,
    p_cutoff = p_cutoff
  )

  .write_suppa_xlsx(
    sig_event_list = sig_event_list,
    sig_gene_list = sig_gene_list,
    summary_tbl = summary_tbl,
    out_xlsx = out_xlsx
  )

  out <- list(
    meta = list(
      project_root = PROJECT_ROOT,
      suppa_dir = suppa_dir,
      ref_level = ref_level,
      timepoints = timepoints,
      p_cutoff = p_cutoff
    ),
    by_tp = by_tp,
    results = list(
      event_full_by_tp = full_list,
      event_sig_by_tp = sig_event_list,
      gene_sig_by_tp = sig_gene_list,
      event_full_all = full_all,
      event_sig_all = sig_all,
      gene_sig_all = sig_genes_all
    ),
    summary = summary_tbl,
    paths = list(rds = out_rds, xlsx = out_xlsx)
  )

  saveRDS(out, out_rds)
  out
}
