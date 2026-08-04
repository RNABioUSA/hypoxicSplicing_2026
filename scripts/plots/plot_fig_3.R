# ============================================================
# plot_fig_3.R
# ============================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
if (!nzchar(PROJECT_ROOT)) stop("PROJECT_ROOT env var not set.")

HELPERS_FILE <- file.path(PROJECT_ROOT, "scripts/utils/helpers.R")
if (!file.exists(HELPERS_FILE)) stop("helpers.R not found at: ", HELPERS_FILE)
source(HELPERS_FILE)

COLORS_FILE <- file.path(PROJECT_ROOT, "scripts/utils/color_palette.R")
if (!file.exists(COLORS_FILE)) stop("color_palette.R not found at: ", COLORS_FILE)
source(COLORS_FILE)

# ---------------------------------------
# Internal Helpers
# ----------------------------------------

.get_condition_labels <- function(conditions) {
  condition_label_map <- c(
    "C1"  = "Normoxia",
    "H1"  = "Hypoxia (1H)",
    "H3"  = "Hypoxia (3H)",
    "H24" = "Hypoxia (24H)"
  )

  labels <- condition_label_map[conditions]

  if (any(is.na(labels))) {
    stop("Unknown condition(s): ", paste(conditions[is.na(labels)], collapse = ", "))
  }

  labels
}

.get_splicing_shades <- function(tool, labels = c("light", "base", "dark")) {
  if (!tool %in% names(splicing_tool_colors)) {
    stop("Unknown splicing tool color: ", tool)
  }

  .generate_shaded_palette(
    base_colors = stats::setNames(unname(splicing_tool_colors[tool]), tool),
    labels = labels
  )[[tool]]
}

.get_highlight_shade <- function() {
  .get_splicing_shades("DRIMSeq")["light"]
}

.get_model_shades <- function() {
  .get_splicing_shades("SUPPA2")
}

.get_gene_locus_title <- function(sashimi_data) {
  paste0(
    sashimi_data$gene_name, " Locus (",
    as.character(GenomicRanges::seqnames(
      sashimi_data$flatExonsByGene[[sashimi_data$returned_gene_key]]
    )[1]),
    ")"
  )
}

.resolve_splicejam_gene <- function(gene, tx2geneDF) {
  if (!is.character(gene) || length(gene) != 1 || is.na(gene) || !nzchar(gene)) {
    stop("gene must be a single non-empty character string.")
  }

  if (is.null(tx2geneDF) || !is.data.frame(tx2geneDF)) {
    stop("tx2geneDF must be a data.frame.")
  }

  req_cols <- c("gene_id", "gene_name")
  missing_cols <- setdiff(req_cols, colnames(tx2geneDF))
  if (length(missing_cols) > 0) {
    stop("tx2geneDF missing required column(s): ", paste(missing_cols, collapse = ", "))
  }

  gene <- as.character(gene)

  # Already an Ensembl gene ID
  if (grepl("^ENSG", gene)) {
    gene_ids <- unique(as.character(tx2geneDF$gene_id[tx2geneDF$gene_id == gene]))
  } else {
    gene_ids <- unique(as.character(tx2geneDF$gene_id[tx2geneDF$gene_name == gene]))
  }

  gene_ids <- gene_ids[!is.na(gene_ids) & nzchar(gene_ids)]

  if (length(gene_ids) == 0) {
    stop("Could not resolve gene to Ensembl ID: ", gene)
  }

  if (length(gene_ids) > 1) {
    stop("Gene maps to multiple Ensembl IDs: ", gene, " -> ", paste(gene_ids, collapse = ", "))
  }

  gene_ids[[1]]
}

.require_two_conditions <- function(conditions, caller = deparse(sys.call(-1)[[1]])) {
  if (!is.character(conditions) || length(conditions) != 2) {
    stop(caller, " requires exactly two conditions, e.g. c('C1', 'H3').")
  }
}

.make_sashimi_files_df <- function(condition, strand) {
  star_dir <- file.path(.get_results_dir(), "counts", "star")

  junc_file <- file.path(star_dir, condition, paste0(condition, ".SJ.cpm.out.tab"))
  bw_file <- file.path(star_dir, condition, paste0(condition, ".bw"))

  if (!file.exists(junc_file)) stop("Missing junction file: ", junc_file)
  if (!file.exists(bw_file)) stop("Missing bigWig file: ", bw_file)

  filesDF <- data.frame(
    stringsAsFactors = FALSE,
    check.names = FALSE,
    url = c(junc_file, bw_file),
    type = c("junction", "bw"),
    sample_id = c(condition, condition),
    scale_factor = 1
  )

  if (strand == "-") {
    filesDF$scale_factor[filesDF$type == "bw"] <- -1
  }

  filesDF
}

.combine_filesDF <- function(sashimi_data, conditions) {
  missing <- setdiff(conditions, names(sashimi_data$filesDF_by_condition))
  if (length(missing) > 0) {
    stop("No filesDF found for condition(s): ", paste(missing, collapse = ", "))
  }

  do.call(rbind, sashimi_data$filesDF_by_condition[conditions])
}

.add_zoom_highlight <- function(p, zoom_region) {
  if (is.null(zoom_region)) {
    return(p)
  }

  highlight_col <- unname(.get_highlight_shade())

  p +
    ggplot2::annotate(
      "rect",
      xmin = zoom_region[1],
      xmax = zoom_region[2],
      ymin = -Inf,
      ymax = Inf,
      fill = highlight_col,
      alpha = 0.35,
      color = highlight_col,
      linewidth = 0.5
    )
}

.rescale <- function(x) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))

  if (sum(ok) == 0) {
    return(out)
  }
  if (length(unique(x[ok])) == 1) {
    out[ok] <- 1
    return(out)
  }

  rng <- range(x[ok], na.rm = TRUE)
  out[ok] <- (x[ok] - rng[1]) / (rng[2] - rng[1])
  out
}

.first_present <- function(tbl, candidates) {
  x <- candidates[candidates %in% colnames(tbl)]
  if (length(x) == 0) {
    return(NULL)
  }
  x[1]
}

.summarize_gene_metric <- function(tbl, gene_col = "ensgene", sig_col = NULL, effect_col = NULL) {
  if (is.null(tbl) || !is.data.frame(tbl) || nrow(tbl) == 0) {
    return(data.frame(
      ensgene = character(),
      best_sig = numeric(),
      best_effect = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  if (!(gene_col %in% colnames(tbl))) {
    stop("Expected gene column '", gene_col, "' not found.")
  }

  keep_cols <- unique(c(gene_col, sig_col, effect_col))
  keep_cols <- keep_cols[!is.na(keep_cols) & !is.null(keep_cols) & nzchar(keep_cols)]
  keep_cols <- keep_cols[keep_cols %in% colnames(tbl)]

  x <- tbl[, keep_cols, drop = FALSE]
  x <- x[!is.na(x[[gene_col]]) & x[[gene_col]] != "", , drop = FALSE]

  if (nrow(x) == 0) {
    return(data.frame(
      ensgene = character(),
      best_sig = numeric(),
      best_effect = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  x[[gene_col]] <- as.character(x[[gene_col]])

  if (!is.null(sig_col) && sig_col %in% colnames(x)) {
    x[[sig_col]] <- suppressWarnings(as.numeric(x[[sig_col]]))
  }
  if (!is.null(effect_col) && effect_col %in% colnames(x)) {
    x[[effect_col]] <- suppressWarnings(as.numeric(x[[effect_col]]))
  }

  split_idx <- split(seq_len(nrow(x)), x[[gene_col]])

  best_sig <- if (!is.null(sig_col) && sig_col %in% colnames(x)) {
    vapply(split_idx, function(idx) {
      vals <- x[[sig_col]][idx]
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) NA_real_ else min(vals)
    }, numeric(1))
  } else {
    setNames(rep(NA_real_, length(split_idx)), names(split_idx))
  }

  best_effect <- if (!is.null(effect_col) && effect_col %in% colnames(x)) {
    vapply(split_idx, function(idx) {
      vals <- x[[effect_col]][idx]
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) NA_real_ else max(abs(vals))
    }, numeric(1))
  } else {
    setNames(rep(NA_real_, length(split_idx)), names(split_idx))
  }

  out <- data.frame(
    ensgene = names(split_idx),
    best_sig = as.numeric(best_sig),
    best_effect = as.numeric(best_effect),
    stringsAsFactors = FALSE
  )

  rownames(out) <- NULL
  out
}

.make_tx2gene_df <- function(txdb_info) {
  .check_pkg("splicejam")

  splicejam::makeTx2geneFromGtf(
    txdb_info$gtf_path,
    geneAttrNames = c("gene_id", "gene_name", "gene_type"),
    txAttrNames = c("transcript_id", "transcript_type"),
    geneFeatureType = "gene",
    txFeatureType = c("transcript", "mRNA")
  )
}

.prep_sashimi_data <- function(
  gene,
  conditions = c("C1", "H3"),
  txdb_info = NULL,
  tx2geneDF = NULL,
  se = NULL,
  force = FALSE
) {
  .check_pkg(c(
    "splicejam",
    "GenomicFeatures",
    "GenomicRanges",
    "S4Vectors"
  ))

  if (is.null(txdb_info)) {
    paths <- .get_runtime_paths()
    env_file <- paths$env_file
    resources_dir <- paths$resources_dir

    txdb_info <- .get_txdb(
      env_file = env_file,
      resources_dir = resources_dir,
      overwrite = FALSE
    )
  }

  txdb <- txdb_info$txdb

  if (is.null(tx2geneDF)) {
    tx2geneDF <- .make_tx2gene_df(txdb_info)
  }

  gene_id <- unique(as.character(tx2geneDF$gene_id[tx2geneDF$gene_name == gene]))
  gene_id <- gene_id[nzchar(gene_id)]

  if (length(gene_id) == 0) {
    stop("Could not resolve gene symbol to Ensembl ID: ", gene)
  }
  if (length(gene_id) > 1) {
    stop("Gene symbol maps to multiple Ensembl IDs: ", gene)
  }

  det <- .get_detected_transcripts_for_gene(
    gene = gene_id,
    conditions = conditions,
    tx2geneDF = tx2geneDF,
    se = se,
    force = force
  )

  gene_symbol <- det$gene_name
  gene_name <- det$gene_name
  detectedTx <- det$detectedTx

  if (length(detectedTx) == 0) {
    stop("No detected transcripts passed filtering for gene: ", gene)
  }

  exonsByTx <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
  S4Vectors::values(exonsByTx@unlistData)$feature_type <- "exon"
  S4Vectors::values(exonsByTx@unlistData)$subclass <- "exon"

  cdsByTx <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)
  S4Vectors::values(cdsByTx@unlistData)$feature_type <- "cds"
  S4Vectors::values(cdsByTx@unlistData)$subclass <- "cds"

  flatExonsByGene <- splicejam::flattenExonsBy(
    exonsByTx = exonsByTx,
    cdsByTx = cdsByTx,
    by = "gene",
    genes = gene_symbol,
    tx2geneDF = tx2geneDF,
    geneColname = "gene_name",
    txColname = "transcript_id",
    verbose = FALSE
  )

  if (length(flatExonsByGene) == 0) {
    stop("flattenExonsBy(by='gene') returned no models for gene: ", gene_symbol)
  }

  returned_gene_key <- names(flatExonsByGene)[1]

  flatExonsByTx <- splicejam::flattenExonsBy(
    exonsByTx = exonsByTx,
    cdsByTx = cdsByTx,
    by = "tx",
    genes = gene_symbol,
    tx2geneDF = tx2geneDF,
    geneColname = "gene_name",
    txColname = "transcript_id",
    detectedTx = detectedTx
  )

  strand <- as.character(
    GenomicRanges::strand(flatExonsByGene[[returned_gene_key]])[1]
  )

  filesDF_by_condition <- setNames(
    lapply(conditions, .make_sashimi_files_df, strand = strand),
    conditions
  )

  list(
    gene = gene,
    gene_id = gene_id,
    gene_symbol = gene_symbol,
    gene_name = gene_name,
    detectedTx = detectedTx,
    tx2geneDF = tx2geneDF,
    flatExonsByGene = flatExonsByGene,
    flatExonsByTx = flatExonsByTx,
    returned_gene_key = returned_gene_key,
    condition_color_map = .get_condition_shaded_colors(conditions, shade = "base"),
    conditions = conditions,
    filesDF_by_condition = filesDF_by_condition
  )
}

.parse_suppa_event_id <- function(event_id) {
  if (!is.character(event_id) || length(event_id) != 1 || is.na(event_id) || !nzchar(event_id)) {
    stop("event_id must be a single non-empty string.")
  }

  event_id <- trimws(gsub("\r", "", gsub('"', "", event_id, fixed = TRUE), fixed = TRUE))

  # Split "ENSG...;TYPE:chr...:<payload>:<strand>"
  parts <- strsplit(event_id, ";", fixed = TRUE)[[1]]
  if (length(parts) != 2) {
    stop("Could not parse event_id: ", event_id)
  }

  gene_part <- parts[1]
  rest <- parts[2]

  rest_parts <- strsplit(rest, ":", fixed = TRUE)[[1]]
  if (length(rest_parts) < 3) {
    stop("Could not parse event_id payload: ", event_id)
  }

  event_type <- rest_parts[1]
  chr <- rest_parts[2]
  strand <- rest_parts[length(rest_parts)]
  payload <- rest_parts[3:(length(rest_parts) - 1)]

  list(
    gene = gene_part,
    event_type = event_type,
    chr = chr,
    strand = strand,
    payload = payload,
    event_id = event_id
  )
}

.extract_event_zoom_coords <- function(event_id, pad = 250) {
  ev <- .parse_suppa_event_id(event_id)
  type <- ev$event_type
  p <- ev$payload

  .parse_coord_token <- function(x) {
    vals <- as.numeric(strsplit(x, "-", fixed = TRUE)[[1]])
    vals[is.finite(vals)]
  }

  .token_start <- function(x) {
    .parse_coord_token(x)[1]
  }

  .token_end <- function(x) {
    vals <- .parse_coord_token(x)
    vals[length(vals)]
  }

  # Returns the local event window BEFORE padding
  core <- switch(type,

    # SE: e1-s2:e2-s3  -> zoom around skipped exon body (s2 to e2) plus downstream splice anchor s3 if desired
    SE = {
      if (length(p) != 2) stop("Unexpected SE payload in: ", event_id)
      c(.token_end(p[1]), .token_end(p[2])) # s2 to s3
    },

    # MX: e1-s2:e2-s4:e1-s3:e3-s4 -> local alternative exon region between s2 and e3
    MX = {
      if (length(p) != 4) stop("Unexpected MX payload in: ", event_id)
      c(.token_end(p[1]), .token_end(p[4])) # s2 to e3
    },

    # A5: e2-s3:e1-s3 -> zoom around alternative donor exon ends to common acceptor
    A5 = {
      if (length(p) != 2) stop("Unexpected A5 payload in: ", event_id)
      vals <- c(.parse_coord_token(p[1]), .parse_coord_token(p[2]))
      range(vals, na.rm = TRUE)
    },

    # A3: e1-s2:e1-s3 -> zoom around common donor to alternative acceptors
    A3 = {
      if (length(p) != 2) stop("Unexpected A3 payload in: ", event_id)
      vals <- c(.parse_coord_token(p[1]), .parse_coord_token(p[2]))
      range(vals, na.rm = TRUE)
    },

    # RI: s1:e1-s2:e2 -> zoom retained intron body from e1 to s2
    RI = {
      if (length(p) != 3) stop("Unexpected RI payload in: ", event_id)
      c(.token_start(p[2]), .token_end(p[2])) # e1 to s2
    },

    # AF: s1:e1-s3:s2:e2-s3 -> zoom local competing first-exon junction region, s2 to s3
    AF = {
      if (length(p) != 4) stop("Unexpected AF payload in: ", event_id)
      c(.token_start(p[3]), .token_end(p[4])) # s2 to s3
    },

    # AL: e1-s2:e2:e1-s3:e3 -> zoom local competing last-exon junction region, s2 to e3
    AL = {
      if (length(p) != 4) stop("Unexpected AL payload in: ", event_id)
      c(.token_end(p[1]), .token_start(p[4])) # s2 to e3
    },

    # fallback: broad range of all parsed coords
    {
      vals <- unlist(lapply(p, .parse_coord_token))
      vals <- vals[is.finite(vals)]
      if (length(vals) < 2) stop("Could not parse coordinates from: ", event_id)
      range(vals, na.rm = TRUE)
    }
  )

  core <- sort(as.numeric(core))
  c(core[1] - pad, core[2] + pad)
}

.get_zoom_region <- function(
  suppa_results,
  gene,
  tp,
  event_id = NULL,
  pad = 250
) {
  suppa_tbl <- suppa_results$by_tp[[tp]]$sig

  if (is.null(suppa_tbl) || !is.data.frame(suppa_tbl)) {
    stop("suppa_results$by_tp[['", tp, "']]$sig not found or not a data.frame.")
  }

  req_cols <- c("symbol", "event_id", "dpsi")
  missing_cols <- setdiff(req_cols, colnames(suppa_tbl))
  if (length(missing_cols) > 0) {
    stop(
      "SUPPA sig table missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  rows <- suppa_tbl[suppa_tbl$symbol == gene, , drop = FALSE]
  if (nrow(rows) == 0) {
    stop("No SUPPA events found for gene: ", gene, " at timepoint: ", tp)
  }

  rows$event_id <- trimws(gsub("\r", "", gsub('"', "", rows$event_id, fixed = TRUE), fixed = TRUE))

  if (!is.null(event_id)) {
    event_id <- trimws(gsub("\r", "", gsub('"', "", event_id, fixed = TRUE), fixed = TRUE))

    strip_event_gene_version <- function(x) {
      sub("^(ENSG[0-9]+)\\.[0-9]+;", "\\1;", x)
    }

    rows_match <- rows[rows$event_id == event_id, , drop = FALSE]

    if (nrow(rows_match) == 0) {
      rows_match <- rows[
        strip_event_gene_version(rows$event_id) == strip_event_gene_version(event_id), ,
        drop = FALSE
      ]
    }

    if (nrow(rows_match) == 0) {
      stop("Requested event_id not found for gene ", gene, " at timepoint ", tp, ": ", event_id)
    }

    rows <- rows_match
  }

  rows$dpsi_num <- suppressWarnings(as.numeric(rows$dpsi))
  rows <- rows[order(-abs(rows$dpsi_num)), , drop = FALSE]

  best <- rows[1, , drop = FALSE]

  best_event_id <- as.character(best$event_id[[1]])
  best_event_type <- if ("event_type" %in% colnames(best)) {
    as.character(best$event_type[[1]])
  } else {
    .parse_suppa_event_id(best_event_id)$event_type
  }

  zoom_region <- .extract_event_zoom_coords(best_event_id, pad = pad)

  list(
    event_row = best,
    event_id = best_event_id,
    event_type = best_event_type,
    dpsi = as.numeric(best$dpsi_num[[1]]),
    zoom_region = zoom_region
  )
}

.read_suppa_psi <- function(condition) {
  psi_file <- file.path(
    .get_results_dir(),
    "analysis",
    "suppa",
    "standard",
    "C1",
    "events",
    paste0("events_", condition, ".psi")
  )

  if (!file.exists(psi_file)) {
    stop("SUPPA PSI file not found: ", psi_file)
  }

  lines <- readLines(psi_file, warn = FALSE)
  if (length(lines) < 2) {
    stop("SUPPA PSI file appears empty or malformed: ", psi_file)
  }

  # First line contains only sample names
  sample_names <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  sample_names <- trimws(gsub("\r", "", sample_names, fixed = TRUE))

  # Remaining lines contain event_id + PSI values
  x <- utils::read.delim(
    text = paste(lines[-1], collapse = "\n"),
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  expected_cols <- length(sample_names) + 1
  if (ncol(x) != expected_cols) {
    stop(
      "Unexpected number of columns in SUPPA PSI file.\n",
      "File: ", psi_file, "\n",
      "Expected: ", expected_cols, "\n",
      "Observed: ", ncol(x)
    )
  }

  colnames(x) <- c("event_id", sample_names)

  x$event_id <- trimws(gsub("\r", "", x$event_id, fixed = TRUE))
  x
}

.get_detected_transcripts_for_gene <- function(
  gene,
  conditions = c("C1", "H3"),
  tx2geneDF,
  se,
  min_tpm = 1,
  top_n = 6,
  force = FALSE
) {
  gene_id <- .resolve_splicejam_gene(gene, tx2geneDF)

  if (is.null(se)) {
    paths <- .get_runtime_paths()
    coldata_tsv <- paths$coldata_tsv
    env_file <- paths$env_file
    resources_dir <- paths$resources_dir
    indexDir <- paths$indexDir

    se <- .get_se(
      coldata_tsv = coldata_tsv,
      env_file = env_file,
      resources_dir = resources_dir,
      indexDir = indexDir,
      force = force
    )
  }

  gene_tx <- unique(as.character(tx2geneDF$transcript_id[tx2geneDF$gene_id == gene_id]))
  if (length(gene_tx) == 0) {
    stop("No transcripts found in tx2geneDF for gene: ", gene_id)
  }

  se_conditions <- as.character(SummarizedExperiment::colData(se)$condition)
  sample_keep <- se_conditions %in% conditions
  if (!any(sample_keep)) {
    stop("No samples matched requested conditions: ", paste(conditions, collapse = ", "))
  }

  assay_names <- SummarizedExperiment::assayNames(se)
  abundance_name <- assay_names[tolower(assay_names) == "abundance"][1]
  if (is.na(abundance_name) || !nzchar(abundance_name)) {
    stop(
      "Could not find transcript-level abundance assay in se.\nAvailable assays: ",
      paste(assay_names, collapse = ", ")
    )
  }

  abundance_mat <- SummarizedExperiment::assay(se, abundance_name)
  tx_keep <- rownames(abundance_mat) %in% gene_tx
  if (!any(tx_keep)) {
    stop("None of the gene transcripts were found in transcript-level abundance matrix.")
  }

  abundance_sub <- abundance_mat[tx_keep, sample_keep, drop = FALSE]
  sample_conditions <- se_conditions[sample_keep]

  summary_df <- data.frame(
    transcript_id = rownames(abundance_sub),
    mean_tpm = rowMeans(abundance_sub, na.rm = TRUE),
    max_tpm = apply(abundance_sub, 1, max, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  for (cond in unique(sample_conditions)) {
    cond_idx <- sample_conditions == cond
    summary_df[[paste0("mean_tpm_", cond)]] <- rowMeans(abundance_sub[, cond_idx, drop = FALSE], na.rm = TRUE)
    summary_df[[paste0("max_tpm_", cond)]] <- apply(abundance_sub[, cond_idx, drop = FALSE], 1, max, na.rm = TRUE)
  }

  summary_df <- summary_df[order(-summary_df$max_tpm, -summary_df$mean_tpm), , drop = FALSE]

  detected_tx <- summary_df$transcript_id[summary_df$max_tpm >= min_tpm]
  if (!is.null(top_n) && length(detected_tx) > top_n) {
    detected_tx <- detected_tx[seq_len(top_n)]
  }

  list(
    gene_id = gene_id,
    gene_name = unique(as.character(tx2geneDF$gene_name[tx2geneDF$gene_id == gene_id]))[1],
    detectedTx = detected_tx,
    summary = summary_df,
    abundance = abundance_sub
  )
}

.plot_candidate_gene_splicejam_model <- function(
  sashimi_data,
  zoom_region = NULL,
  labelExons = FALSE,
  exonLabelSize = 3,
  save_plot = FALSE,
  plot_width = 12,
  plot_height = 4.5
) {
  .check_pkg(c(
    "splicejam",
    "ggplot2",
    "cowplot"
  ))

  gene <- sashimi_data$gene
  gene_id <- sashimi_data$gene_id
  gene_symbol <- sashimi_data$gene_symbol
  gene_name <- sashimi_data$gene_name
  tx2geneDF <- sashimi_data$tx2geneDF
  flatExonsByGene <- sashimi_data$flatExonsByGene
  flatExonsByTx <- sashimi_data$flatExonsByTx
  returned_gene_key <- sashimi_data$returned_gene_key

  annot_shades <- .get_model_shades()
  highlight_shade <- .get_highlight_shade()

  gg_model <- suppressWarnings(
    splicejam::gene2gg(
      gene = returned_gene_key,
      flatExonsByGene = flatExonsByGene,
      flatExonsByTx = flatExonsByTx,
      tx2geneDF = tx2geneDF,
      labelExons = labelExons,
      exonLabelSize = exonLabelSize
    )
  ) +
    ggplot2::ggtitle(gene_name)

  if (is.data.frame(gg_model$data) && "color_by" %in% names(gg_model$data)) {
    gg_model$data$color_by <- sub(
      "^.*_(noncds|cds|gap)$",
      "\\1",
      as.character(gg_model$data$color_by)
    )
  }

  for (i in seq_along(gg_model$layers)) {
    if (is.data.frame(gg_model$layers[[i]]$data) &&
      nrow(gg_model$layers[[i]]$data) > 0 &&
      "color_by" %in% names(gg_model$layers[[i]]$data)) {
      gg_model$layers[[i]]$data$color_by <- sub(
        "^.*_(noncds|cds|gap)$",
        "\\1",
        as.character(gg_model$layers[[i]]$data$color_by)
      )
    }
  }

  if (length(gg_model$scales$scales) > 0) {
    keep <- vapply(gg_model$scales$scales, function(s) {
      !any(s$aesthetics %in% c("fill", "colour", "color"))
    }, logical(1))
    gg_model$scales$scales <- gg_model$scales$scales[keep]
  }

  if (!is.null(zoom_region)) {
    gg_model <- gg_model +
      ggplot2::annotate(
        "rect",
        xmin = zoom_region[1],
        xmax = zoom_region[2],
        ymin = -Inf,
        ymax = Inf,
        fill = unname(highlight_shade),
        alpha = 0.35,
        color = unname(highlight_shade),
        linewidth = 0.5
      )
  }

  gg_model <- gg_model +
    ggplot2::scale_fill_manual(
      values = c(
        cds = unname(annot_shades["dark"]),
        noncds = unname(annot_shades["light"]),
        gap = unname(annot_shades["dark"])
      ),
      breaks = c("cds", "noncds", "gap")
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        cds = unname(annot_shades["dark"]),
        noncds = unname(annot_shades["light"]),
        gap = unname(annot_shades["dark"])
      ),
      breaks = c("cds", "noncds", "gap")
    )

  pb0 <- ggplot2::ggplot_build(gg_model)
  pp0 <- pb0$layout$panel_params[[1]]
  y_breaks <- pp0$y$breaks
  y_labels <- pp0$y$get_labels()

  y_labels_new <- y_labels
  y_labels_new <- sub(
    pattern = paste0("^", gene_symbol, "_"),
    replacement = paste0(gene_id, "_"),
    x = y_labels_new
  )
  y_labels_new[y_labels_new == gene_symbol] <- gene_id
  y_labels_new[y_labels_new == returned_gene_key] <- gene_id

  gg_model <- gg_model +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels_new
    )

  gene_idx <- which(y_labels_new == gene_id)[1]
  if (!is.na(gene_idx) && length(y_breaks) >= 2) {
    gene_y <- y_breaks[gene_idx]
    other_y <- y_breaks[seq_along(y_breaks) != gene_idx]
    nearest_y <- other_y[which.min(abs(other_y - gene_y))]
    sep_y <- mean(c(gene_y, nearest_y))

    gg_model <- gg_model +
      ggplot2::geom_hline(
        yintercept = sep_y,
        linewidth = 0.75,
        color = "black"
      )
  }

  pb_model <- ggplot2::ggplot_build(gg_model)
  y_vals <- unlist(lapply(pb_model$data, function(d) {
    y_cols <- intersect(c("y", "ymin", "ymax", "yend"), names(d))
    if (length(y_cols) == 0) {
      return(NULL)
    }
    unlist(d[y_cols], use.names = FALSE)
  }))
  y_vals <- y_vals[is.finite(y_vals)]

  if (length(y_vals) > 0) {
    y_rng <- range(y_vals, na.rm = TRUE)
    y_pad <- 0.5

    gg_model <- gg_model +
      ggplot2::coord_cartesian(
        ylim = c(y_rng[1] - y_pad, y_rng[2] + y_pad),
        expand = TRUE
      )
  }

  x_title <- paste0(
    gene_name, " Locus (",
    as.character(GenomicRanges::seqnames(flatExonsByGene[[returned_gene_key]])[1]),
    ")"
  )

  gg_model <- gg_model +
    ggplot2::labs(x = x_title) +
    ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 18),
      plot.margin = ggplot2::margin(5.5, 12, 2, 12)
    )

  p_out <- gg_model

  saved_file <- NULL
  if (isTRUE(save_plot)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_3")
    gene_dir <- file.path(outdir, gene)
    dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

    outfile <- file.path(gene_dir, paste0("fig_3_A_model_", gene, ".pdf"))
    ggplot2::ggsave(outfile, p_out, width = plot_width, height = plot_height)
    saved_file <- outfile
  }

  list(
    plot = p_out,
    base_plot = gg_model,
    saved_file = saved_file
  )
}

.plot_candidate_gene_splicejam_sashimi <- function(
  sashimi_data,
  conditions = c("C1", "H3"),
  zoom_region = NULL,
  minJunctionScore = 1,
  show = c("coverage", "junction", "junctionLabels"),
  show_x_title = TRUE,
  save_plot = FALSE,
  plot_width = 9,
  plot_height = 5 * length(conditions)
) {
  .check_pkg(c("splicejam", "ggplot2", "ggh4x", "jamba"))

  gene <- sashimi_data$gene
  color_sub <- .get_condition_shaded_colors(conditions, shade = "base")
  condition_labels <- .get_condition_labels(conditions)
  filesDF <- .combine_filesDF(sashimi_data, conditions)

  sh <- suppressWarnings(
    splicejam::prepareSashimi(
      gene = sashimi_data$returned_gene_key,
      flatExonsByGene = sashimi_data$flatExonsByGene,
      minJunctionScore = minJunctionScore,
      sample_id = conditions,
      filesDF = filesDF
    )
  )

  gg_sashimi <- suppressMessages(
    suppressWarnings(
      splicejam::plotSashimi(
        sh,
        junc_color = jamba::alpha2col(unname(color_sub), 0.45),
        junc_fill = jamba::alpha2col(unname(color_sub), 0.35),
        show = show,
        fill_scheme = "sample_id",
        color_sub = color_sub,
        apply_facet = TRUE,
        verbose = FALSE
      )
    )
  ) +
    ggplot2::ggtitle(NULL)

  gg_sashimi <- .add_zoom_highlight(gg_sashimi, zoom_region)

  gg_sashimi <- gg_sashimi +
    ggh4x::facet_grid2(
      rows = ggplot2::vars(sample_id),
      scales = "free_y",
      labeller = ggplot2::labeller(
        sample_id = stats::setNames(condition_labels, conditions)
      ),
      strip = ggh4x::strip_themed(
        background_y = ggh4x::elem_list_rect(
          fill = unname(color_sub),
          color = NA
        ),
        text_y = ggh4x::elem_list_text(
          size = 15,
          color = "black"
        )
      )
    )

  x_title <- .get_gene_locus_title(sashimi_data)
  y_title <- "Read Depth (CPM)"

  if (length(gg_sashimi$scales$scales) > 0) {
    for (i in seq_along(gg_sashimi$scales$scales)) {
      aes_i <- gg_sashimi$scales$scales[[i]]$aesthetics
      if ("x" %in% aes_i) {
        gg_sashimi$scales$scales[[i]]$name <- if (isTRUE(show_x_title)) x_title else NULL
      }
      if ("y" %in% aes_i) {
        gg_sashimi$scales$scales[[i]]$name <- y_title
      }
    }
  }

  gg_sashimi <- gg_sashimi +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 15),
      plot.margin = ggplot2::margin(5.5, 12, 5.5, 12),
      axis.title.x = ggplot2::element_text(size = 15, margin = ggplot2::margin(t = 12.5)),
      axis.title.y = ggplot2::element_text(size = 15)
    )

  saved_file <- NULL
  if (isTRUE(save_plot)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_3")
    gene_dir <- file.path(outdir, gene)
    dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

    outfile <- file.path(
      gene_dir,
      paste0("fig_3_A_sashimi_", paste(conditions, collapse = "_"), "_", gene, ".pdf")
    )
    ggplot2::ggsave(outfile, gg_sashimi, width = plot_width, height = plot_height)
    saved_file <- outfile
  }

  list(
    plot = gg_sashimi,
    base_plot = gg_sashimi,
    sh = sh,
    saved_file = saved_file
  )
}

.plot_candidate_gene_splicejam_sashimi_zoom <- function(
  sashimi_res,
  zoom_region,
  event_id = NULL,
  zoom_size = 1,
  save_plot = TRUE,
  plot_width = 9,
  plot_height = 6.5
) {
  .check_pkg(c(
    "ggplot2",
    "ggforce",
    "jamba"
  ))

  sh <- sashimi_res$sh
  p_orig <- sashimi_res$base_plot
  condition <- unique(as.character(sh$df$sample_id))[1]

  base_col <- unname(.get_condition_colors(condition))
  highlight_shade <- .get_highlight_shade()
  zoom_fill <- unname(highlight_shade)
  zoom_xlim_comp <- sh$ref2c$transform(sort(zoom_region))

  pb <- ggplot2::ggplot_build(p_orig)
  cov_df <- pb$data[[1]]
  junc_df <- pb$data[[2]]

  label_df <- NULL
  if (length(pb$data) >= 3) {
    label_df <- pb$data[[3]]
    label_df <- label_df[, intersect(c("x", "y", "label"), names(label_df)), drop = FALSE]

    if (nrow(label_df) == 0 || !all(c("x", "y", "label") %in% names(label_df))) {
      label_df <- NULL
    }
  }

  cov_df2 <- cov_df[, intersect(c("x", "y", "group"), names(cov_df)), drop = FALSE]
  junc_df2 <- junc_df[, intersect(c("x", "y", "group"), names(junc_df)), drop = FALSE]

  p_combo <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = cov_df2,
      ggplot2::aes(x = x, y = y, group = group),
      inherit.aes = FALSE,
      fill = base_col,
      color = NA,
      alpha = 0.9
    ) +
    ggplot2::geom_polygon(
      data = junc_df2,
      ggplot2::aes(x = x, y = y, group = group),
      inherit.aes = FALSE,
      fill = jamba::alpha2col(base_col, 0.35),
      color = jamba::alpha2col(base_col, 0.55)
    )

  if (!is.null(label_df)) {
    p_combo <- p_combo +
      ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        size = 3
      )
  }

  p_combo <- p_combo +
    ggforce::facet_zoom(
      xlim = zoom_xlim_comp,
      horizontal = FALSE,
      zoom.size = zoom_size,
      show.area = TRUE
    ) +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::theme(
      zoom.x = ggplot2::element_rect(
        fill = jamba::alpha2col(zoom_fill, 0.35),
        colour = jamba::alpha2col(zoom_fill, 0.65),
        linewidth = 0.7,
        linetype = 1
      )
    )

  # Recover gene/locus label from original plot
  gene_name <- as.character(sh$gene)
  gene_seqname <- NULL

  x_scales <- Filter(function(s) "x" %in% s$aesthetics, p_orig$scales$scales)
  if (length(x_scales) > 0 && !is.null(x_scales[[1]]$name)) {
    x_title_existing <- x_scales[[1]]$name
    gene_name <- sub(" Locus \\(.*$", "", x_title_existing)

    gene_seqname <- sub("^.* Locus \\((.*)\\)$", "\\1", x_title_existing)
    if (identical(gene_seqname, x_title_existing)) {
      gene_seqname <- NULL
    }
  }

  # Build subtitle from event type / chr + zoom coordinates
  subtitle_txt <- NULL

  if (!is.null(event_id) && nzchar(event_id)) {
    event_short <- sub("^ENSG[0-9.]+;", "", event_id)

    # Extract event type (e.g., A5)
    event_type <- sub("^([A-Z0-9]+):.*$", "\\1", event_short)

    # Extract chromosome (e.g., chr8)
    chr <- sub("^[A-Z0-9]+:(chr[^:]+):.*$", "\\1", event_short)

    # Build clean subtitle using zoom region
    subtitle_txt <- paste0(
      chr, ":",
      format(zoom_region[1], big.mark = ","),
      "-",
      format(zoom_region[2], big.mark = ","),
      " (", event_type, ")"
    )
  } else if (!is.null(zoom_region)) {
    subtitle_txt <- paste0(
      "Zoom: ",
      format(zoom_region[1], big.mark = ","),
      "-",
      format(zoom_region[2], big.mark = ",")
    )
  }

  # Adjust axis labels
  if (is.null(gene_seqname)) {
    x_title <- paste0(gene_name, " Locus")
  } else {
    x_title <- paste0(gene_name, " Locus (", gene_seqname, ")")
  }

  p_combo <- p_combo +
    ggplot2::labs(
      title = gene_name,
      subtitle = subtitle_txt,
      x = x_title,
      y = "Read Depth (CPM)"
    ) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(
        size = 15,
        margin = ggplot2::margin(t = 12.5)
      ),
      axis.title.y = ggplot2::element_text(size = 15),
      plot.title = ggplot2::element_text(size = 18),
      plot.subtitle = ggplot2::element_text(size = 12)
    )

  saved_file <- NULL
  if (isTRUE(save_plot)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_3")
    gene_dir <- file.path(outdir, as.character(sh$gene))
    dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

    outfile <- file.path(
      gene_dir,
      paste0("fig_3_A_sashimi_zoom_", condition, "_", as.character(sh$gene), ".pdf")
    )
    ggplot2::ggsave(outfile, p_combo, width = plot_width, height = plot_height)
    saved_file <- outfile
  }

  list(
    plot = p_combo,
    base_plot = p_combo,
    zoom_xlim_comp = zoom_xlim_comp,
    saved_file = saved_file
  )
}

.plot_candidate_gene_splicejam_sashimi_zoom_only <- function(
  sashimi_data,
  tp = "H3",
  conditions = c("C1", tp),
  zoom_region,
  minJunctionScore = 1,
  show = c("coverage", "junction", "junctionLabels"),
  save_plot = FALSE,
  plot_width = 12,
  plot_height = 6
) {
  .check_pkg(c(
    "splicejam",
    "ggplot2",
    "ggh4x",
    "jamba"
  ))

  .require_two_conditions(conditions, ".plot_candidate_gene_splicejam_sashimi_zoom_only")

  if (!is.numeric(zoom_region) || length(zoom_region) != 2 || any(!is.finite(zoom_region))) {
    stop("zoom_region must be numeric c(start, end) in genomic coordinates.")
  }
  zoom_region <- sort(as.numeric(zoom_region))

  gene <- sashimi_data$gene
  gene_name <- sashimi_data$gene_name
  flatExonsByGene <- sashimi_data$flatExonsByGene
  returned_gene_key <- sashimi_data$returned_gene_key

  color_sub <- .get_condition_shaded_colors(conditions, shade = "base")
  condition_labels <- .get_condition_labels(conditions)
  filesDF <- .combine_filesDF(sashimi_data, conditions)
  x_title <- .get_gene_locus_title(sashimi_data)
  y_title <- "Read Depth (CPM)"

  sh <- suppressWarnings(
    splicejam::prepareSashimi(
      gene = returned_gene_key,
      flatExonsByGene = flatExonsByGene,
      minJunctionScore = minJunctionScore,
      sample_id = conditions,
      filesDF = filesDF
    )
  )

  p_zoom <- suppressMessages(
    suppressWarnings(
      splicejam::plotSashimi(
        sh,
        junc_color = jamba::alpha2col(unname(color_sub), 0.45),
        junc_fill = jamba::alpha2col(unname(color_sub), 0.35),
        show = show,
        fill_scheme = "sample_id",
        label_coords = zoom_region,
        color_sub = color_sub,
        apply_facet = TRUE,
        verbose = FALSE
      )
    )
  ) +
    ggplot2::ggtitle(NULL) +
    ggplot2::coord_cartesian(
      xlim = zoom_region,
      expand = FALSE
    )

  genomic_breaks <- pretty(zoom_region, n = 5)
  genomic_breaks <- genomic_breaks[
    genomic_breaks >= zoom_region[1] & genomic_breaks <= zoom_region[2]
  ]

  p_zoom <- p_zoom +
    ggplot2::scale_x_continuous(
      breaks = genomic_breaks,
      labels = format(genomic_breaks, big.mark = ",", scientific = FALSE)
    )

  p_zoom <- p_zoom +
    ggh4x::facet_grid2(
      rows = ggplot2::vars(sample_id),
      scales = "free_y",
      labeller = ggplot2::labeller(sample_id = stats::setNames(condition_labels, conditions)),
      strip = ggh4x::strip_themed(
        background_y = ggh4x::elem_list_rect(
          fill = unname(color_sub),
          color = NA
        ),
        text_y = ggh4x::elem_list_text(
          size = 15,
          color = "black"
        )
      )
    ) +
    ggplot2::theme_bw(base_size = 15)

  if (length(p_zoom$scales$scales) > 0) {
    for (i in seq_along(p_zoom$scales$scales)) {
      aes_i <- p_zoom$scales$scales[[i]]$aesthetics
      if ("x" %in% aes_i) p_zoom$scales$scales[[i]]$name <- x_title
      if ("y" %in% aes_i) p_zoom$scales$scales[[i]]$name <- y_title
    }
  }

  p_zoom <- p_zoom +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 15),
      plot.margin = ggplot2::margin(5.5, 12, 5.5, 12),
      axis.title.x = ggplot2::element_text(
        size = 15,
        margin = ggplot2::margin(t = 12.5)
      ),
      axis.title.y = ggplot2::element_text(size = 15)
    )

  saved_file <- NULL
  if (isTRUE(save_plot)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_3")
    gene_dir <- file.path(outdir, gene)
    dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

    outfile <- file.path(
      gene_dir,
      paste0("fig_3_A_sashimi_zoom_only_", paste(conditions, collapse = "_"), "_", gene, ".pdf")
    )
    ggplot2::ggsave(outfile, p_zoom, width = plot_width, height = plot_height)
    saved_file <- outfile
  }

  list(
    plot = p_zoom,
    base_plot = p_zoom,
    sh = sh,
    saved_file = saved_file
  )
}

.plot_candidate_gene_splicejam_model_and_sashimi <- function(
  sashimi_data,
  tp = "H3",
  conditions = c("C1", tp),
  zoom_region = NULL,
  minJunctionScore = 1,
  show = c("coverage", "junction", "junctionLabels"),
  rel_heights = c(1, 2),
  show_x_title = FALSE,
  save_plot = FALSE,
  plot_width = 12,
  plot_height = 12
) {
  .check_pkg(c(
    "splicejam",
    "ggplot2",
    "ggh4x",
    "jamba",
    "cowplot"
  ))

  .require_two_conditions(conditions, ".plot_candidate_gene_splicejam_model_and_sashimi")

  gene <- sashimi_data$gene
  gene_name <- sashimi_data$gene_name
  flatExonsByGene <- sashimi_data$flatExonsByGene
  returned_gene_key <- sashimi_data$returned_gene_key

  color_sub <- .get_condition_shaded_colors(conditions, shade = "base")
  highlight_shade <- .get_highlight_shade()

  # Model plot
  p_model <- .plot_candidate_gene_splicejam_model(
    sashimi_data = sashimi_data,
    zoom_region = zoom_region,
    save_plot = FALSE
  )

  # Build one combined two-condition sashimi
  p_sashimi <- .plot_candidate_gene_splicejam_sashimi(
    sashimi_data = sashimi_data,
    conditions = conditions,
    zoom_region = zoom_region,
    minJunctionScore = minJunctionScore,
    show = show,
    show_x_title = show_x_title,
    save_plot = FALSE
  )

  # Stack model + combined sashimi
  p_full <- cowplot::plot_grid(
    p_model$base_plot +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      ) +
      ggplot2::xlab(NULL),
    p_sashimi$base_plot,
    ncol = 1,
    align = "v",
    axis = "lr",
    rel_heights = rel_heights
  )

  saved_file <- NULL
  if (isTRUE(save_plot)) {
    outdir <- file.path(.get_results_dir(), "plots", "fig_3")
    gene_dir <- file.path(outdir, gene)
    dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

    outfile <- file.path(
      gene_dir,
      paste0("fig_3_A_splicejam_full_", gene, ".pdf")
    )
    ggplot2::ggsave(outfile, p_full, width = plot_width, height = plot_height)
    saved_file <- outfile
  }

  list(
    plot = p_full,
    model = p_model,
    sashimi = p_sashimi,
    base_plot = p_sashimi$base_plot,
    sh = p_sashimi$sh,
    saved_file = saved_file
  )
}

.plot_drimseq_ribbon <- function(prop_fit, main = NULL, feature_colors) {
  .check_pkg("ggplot2")

  # Expect prop_fit to contain:
  # feature_id, group, proportion
  stopifnot(all(c("feature_id", "group", "proportion") %in% colnames(prop_fit)))

  prop_fit <- prop_fit[, c("feature_id", "group", "proportion"), drop = FALSE]
  prop_fit$feature_id <- factor(prop_fit$feature_id, levels = levels(prop_fit$feature_id))
  prop_fit$group <- factor(prop_fit$group, levels = levels(prop_fit$group))

  width <- 0.5
  breaks <- levels(prop_fit$feature_id)

  # compute stacked ymin/ymax explicitly, within each group
  stacked_list <- lapply(levels(prop_fit$group), function(g) {
    d <- prop_fit[prop_fit$group == g, , drop = FALSE]

    # enforce identical stacking order in every group
    stack_order <- rev(levels(prop_fit$feature_id))
    d <- d[match(stack_order, d$feature_id), , drop = FALSE]

    d$proportion[is.na(d$proportion)] <- 0
    d$ymax <- cumsum(d$proportion)
    d$ymin <- d$ymax - d$proportion
    d
  })

  stacked_df <- do.call(rbind, stacked_list)
  rownames(stacked_df) <- NULL

  # bar plot
  p <- ggplot2::ggplot(
    stacked_df,
    ggplot2::aes(x = group, y = proportion, fill = feature_id)
  ) +
    ggplot2::geom_col(width = width)

  # ribbons between adjacent groups
  group_levels <- levels(stacked_df$group)

  if (length(group_levels) >= 2) {
    for (i in seq_len(length(group_levels) - 1)) {
      g1 <- group_levels[i]
      g2 <- group_levels[i + 1]

      d1 <- stacked_df[stacked_df$group == g1, , drop = FALSE]
      d2 <- stacked_df[stacked_df$group == g2, , drop = FALSE]

      # align by transcript
      stack_order <- rev(levels(prop_fit$feature_id))
      d1 <- d1[match(stack_order, d1$feature_id), , drop = FALSE]
      d2 <- d2[match(stack_order, d2$feature_id), , drop = FALSE]

      ribbon_df <- do.call(
        rbind,
        lapply(seq_len(nrow(d1)), function(j) {
          data.frame(
            feature_id = as.character(d1$feature_id[j]),
            x = c(
              i + width / 2,
              i + width / 2,
              i + 1 - width / 2,
              i + 1 - width / 2
            ),
            y = c(
              d1$ymin[j],
              d1$ymax[j],
              d2$ymax[j],
              d2$ymin[j]
            ),
            stringsAsFactors = FALSE
          )
        })
      )

      p <- p +
        ggplot2::geom_polygon(
          data = ribbon_df,
          ggplot2::aes(
            x = x,
            y = y,
            fill = feature_id,
            group = feature_id
          ),
          inherit.aes = FALSE,
          alpha = 0.55,
          colour = NA
        )
    }
  }

  p +
    ggplot2::coord_cartesian(ylim = c(-0.05, 1.05)) +
    ggplot2::labs(
      title = main,
      x = "",
      y = "Estimated Transcript Proportion"
    ) +
    ggplot2::theme_bw(base_size = 15)
}

.plot_suppa_psi <- function(
  suppa_results,
  gene,
  tp,
  conditions,
  event_id = NULL,
  zoom_region = NULL,
  show_title = TRUE,
  save_plot = FALSE,
  plot_width = 5,
  plot_height = 5
) {
  .check_pkg("ggplot2")

  paths <- .get_runtime_paths()
  outdir <- paths$outdir

  condition_labels <- .get_condition_labels(conditions)

  if (any(is.na(condition_labels))) {
    stop("Unknown condition(s): ", paste(conditions[is.na(condition_labels)], collapse = ", "))
  }

  if (!tp %in% names(suppa_results$by_tp)) {
    stop("Timepoint not found in suppa_results$by_tp: ", tp)
  }

  suppa_tp <- suppa_results$by_tp[[tp]]
  if (is.null(suppa_tp) || !is.list(suppa_tp) || is.null(suppa_tp$sig)) {
    stop("SUPPA2 significant results table not found for timepoint: ", tp)
  }

  suppa_tbl <- suppa_tp$sig

  req_cols <- c("event_id", "symbol")
  if (!is.data.frame(suppa_tbl) || !all(req_cols %in% colnames(suppa_tbl))) {
    stop("SUPPA2 sig table must contain columns: ", paste(req_cols, collapse = ", "))
  }

  gene_rows <- suppa_tbl[suppa_tbl$symbol == gene, , drop = FALSE]

  if (nrow(gene_rows) == 0) {
    stop("No significant SUPPA2 events found for gene: ", gene, " at timepoint: ", tp)
  }

  gene_rows$event_id <- trimws(gsub("\r", "", gene_rows$event_id, fixed = TRUE))
  gene_rows <- gene_rows[!is.na(gene_rows$event_id) & nzchar(gene_rows$event_id), , drop = FALSE]

  if (nrow(gene_rows) == 0) {
    stop("No valid event_id values found for gene: ", gene, " at timepoint: ", tp)
  }

  if (is.null(event_id)) {
    if ("dpsi" %in% colnames(gene_rows)) {
      gene_rows$dpsi_num <- suppressWarnings(as.numeric(gene_rows$dpsi))
      gene_rows <- gene_rows[order(-abs(gene_rows$dpsi_num)), , drop = FALSE]
    }
    event_id <- gene_rows$event_id[1]
  } else {
    event_id <- trimws(gsub("\r", "", gsub('"', "", event_id, fixed = TRUE), fixed = TRUE))
  }

  subtitle_txt <- NULL

  if (!is.null(event_id) && nzchar(event_id)) {
    event_short <- sub("^ENSG[0-9.]+;", "", event_id)

    event_type <- sub("^([A-Z0-9]+):.*$", "\\1", event_short)
    chr <- sub("^[A-Z0-9]+:(chr[^:]+):.*$", "\\1", event_short)

    if (!is.null(zoom_region) && length(zoom_region) == 2) {
      subtitle_txt <- paste0(
        "SUPPA2 Event\n",
        chr, ":",
        format(zoom_region[1], big.mark = ","),
        "-",
        format(zoom_region[2], big.mark = ","),
        " (", event_type, ")"
      )
    } else {
      subtitle_txt <- event_short
    }
  } else if (!is.null(zoom_region) && length(zoom_region) == 2) {
    subtitle_txt <- paste0(
      "Zoom: ",
      format(zoom_region[1], big.mark = ","),
      "-",
      format(zoom_region[2], big.mark = ",")
    )
  }

  get_event <- function(df, cond, event_id) {
    event_id_norm <- trimws(gsub("\r", "", gsub('"', "", event_id, fixed = TRUE), fixed = TRUE))

    hit <- df[df$event_id == event_id_norm, , drop = FALSE]

    if (nrow(hit) == 0) {
      strip_event_gene_version <- function(x) {
        sub("^(ENSG[0-9]+)\\.[0-9]+;", "\\1;", x)
      }

      hit <- df[
        strip_event_gene_version(df$event_id) == strip_event_gene_version(event_id_norm), ,
        drop = FALSE
      ]
    }

    if (nrow(hit) == 0) {
      cat("\nFirst few PSI event IDs from", cond, ":\n")
      print(utils::head(df$event_id, 5))
      stop("Event not found in ", cond, ": ", event_id)
    }

    if (nrow(hit) > 1) {
      warning("Multiple rows found for event in ", cond, ": ", event_id, ". Using first row.")
      hit <- hit[1, , drop = FALSE]
    }

    hit[, setdiff(colnames(hit), "event_id"), drop = FALSE]
  }

  psi_list <- lapply(conditions, .read_suppa_psi)
  names(psi_list) <- conditions

  psi_event <- mapply(
    get_event,
    df = psi_list,
    cond = conditions,
    MoreArgs = list(event_id = event_id),
    SIMPLIFY = FALSE
  )

  psi_long_list <- lapply(seq_along(conditions), function(i) {
    cond <- conditions[i]
    df <- psi_event[[i]]

    data.frame(
      sample_id = colnames(df),
      PSI = as.numeric(df[1, , drop = TRUE]),
      condition = cond,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })

  psi_long <- do.call(rbind, psi_long_list)
  psi_long$PSI <- suppressWarnings(as.numeric(as.character(psi_long$PSI)))
  psi_long <- psi_long[is.finite(psi_long$PSI), , drop = FALSE]

  if (nrow(psi_long) == 0) {
    stop("No finite PSI values found for event: ", event_id)
  }

  condition_map <- stats::setNames(condition_labels, conditions)
  psi_long$condition <- factor(
    condition_map[psi_long$condition],
    levels = condition_labels
  )

  fill_vals <- .get_condition_shaded_colors(conditions, shade = "base")
  names(fill_vals) <- .get_condition_labels(conditions)

  p <- ggplot2::ggplot(
    psi_long,
    ggplot2::aes(x = condition, y = PSI, fill = condition)
  ) +
    ggplot2::geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.8) +
    ggplot2::geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    ggplot2::scale_fill_manual(values = fill_vals, guide = "none") +
    ggplot2::labs(
      title = if (isTRUE(show_title)) gene else NULL,
      subtitle = subtitle_txt,
      x = "",
      y = "PSI"
    ) +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 18),
      plot.subtitle = ggplot2::element_text(size = 12.5),
      axis.title.y = ggplot2::element_text(size = 15),
      axis.text.x = ggplot2::element_text(size = 15, margin = ggplot2::margin(t = 12.5)),
      axis.text.y = ggplot2::element_text(size = 12)
    )

  saved_file <- NULL
  if (isTRUE(save_plot)) {
    gene_dir <- file.path(outdir, gene)
    dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

    outfile <- file.path(
      gene_dir,
      paste0("fig_3_A_candidate_gene_psi_", gene, ".pdf")
    )

    ggplot2::ggsave(
      outfile,
      p,
      width = plot_width,
      height = plot_height
    )

    saved_file <- outfile
  }

  list(
    plot = p,
    event_id = event_id,
    saved_file = saved_file
  )
}

# --------------------------------------------------
# Candidate Gene Analysis
# --------------------------------------------------

run_splicing_candidate_gene_analysis <- function(
  deseq_results,
  dexseq_results,
  drimseq_results,
  suppa_results,
  timepoints = c("H1", "H3", "H24"),
  focus_mode = c("splicing_only", "shared_with_deg"),
  focus_timepoint = NULL,
  annot_df = NULL
) {
  .check_pkg("dplyr")

  focus_mode <- match.arg(focus_mode)

  paths <- .get_runtime_paths()
  outdir <- paths$outdir

  if (is.null(annot_df)) {
    annot_df <- .get_annot()
  }

  # significant gene sets by timepoint
  message("Starting candidate gene analysis")
  message("focus_mode = ", focus_mode)
  message("timepoints = ", paste(timepoints, collapse = ", "))

  message("Extracting significant gene sets by timepoint")

  deg_genes_by_tp <- .subset_gene_lists_by_tp(
    .get_sig_genes_by_tp(deseq_results$results$gene_sig_by_tp),
    timepoints = timepoints
  )

  dex_genes_by_tp <- .subset_gene_lists_by_tp(
    .get_sig_genes_by_tp(dexseq_results$results$gene_sig_by_tp),
    timepoints = timepoints
  )

  drim_genes_by_tp <- .subset_gene_lists_by_tp(
    .get_sig_genes_by_tp(drimseq_results$results$gene_sig_by_tp),
    timepoints = timepoints
  )

  suppa_genes_by_tp <- .subset_gene_lists_by_tp(
    .get_sig_genes_by_tp(suppa_results$results$gene_sig_by_tp),
    timepoints = timepoints
  )

  message("Building overlap gene sets")

  common_splicing_genes_by_tp <- .intersect_gene_lists_by_tp(
    .intersect_gene_lists_by_tp(dex_genes_by_tp, drim_genes_by_tp),
    suppa_genes_by_tp
  )
  common_splicing_genes_by_tp <- .subset_gene_lists_by_tp(common_splicing_genes_by_tp, timepoints)

  splicing_only_common_genes_by_tp <- .subtract_gene_lists_by_tp(
    common_splicing_genes_by_tp,
    deg_genes_by_tp
  )
  splicing_only_common_genes_by_tp <- .subset_gene_lists_by_tp(splicing_only_common_genes_by_tp, timepoints)

  common_with_deg_genes_by_tp <- .intersect_gene_lists_by_tp(
    common_splicing_genes_by_tp,
    deg_genes_by_tp
  )
  common_with_deg_genes_by_tp <- .subset_gene_lists_by_tp(common_with_deg_genes_by_tp, timepoints)

  candidate_gene_sets_by_tp <- switch(focus_mode,
    splicing_only = splicing_only_common_genes_by_tp,
    shared_with_deg = common_with_deg_genes_by_tp
  )

  message("Done building overlap gene sets")

  for (tp in timepoints) {
    message(
      "", tp,
      ": common_all_3 = ", length(common_splicing_genes_by_tp[[tp]]),
      ", splicing_only = ", length(splicing_only_common_genes_by_tp[[tp]]),
      ", shared_with_deg = ", length(common_with_deg_genes_by_tp[[tp]])
    )
  }

  # candidate rankings by timepoint
  candidate_rankings_by_tp <- setNames(
    lapply(timepoints, function(tp) {
      message("Ranking genes for ", tp)
      candidate_genes <- candidate_gene_sets_by_tp[[tp]]

      message("", tp, ": candidate genes = ", length(candidate_genes))

      if (length(candidate_genes) == 0) {
        return(data.frame(
          timepoint = character(),
          ensgene = character(),
          symbol = character(),
          description = character(),
          dex_best_padj = numeric(),
          dex_best_effect = numeric(),
          drim_best_padj = numeric(),
          drim_best_effect = numeric(),
          suppa_best_pvalue = numeric(),
          suppa_best_dpsi = numeric(),
          suppa_event_types = character(),
          mean_sig_score = numeric(),
          mean_effect_score = numeric(),
          candidate_score = numeric(),
          stringsAsFactors = FALSE
        ))
      }

      # DEXSeq: exon-level table is best for ranking candidate genes
      message("", tp, ": summarizing DEXSeq exon-level metrics")
      dex_time <- system.time({
        dex_tbl <- dexseq_results$results$exon_full_by_tp[[tp]]
        dex_effect_candidates <- c(
          paste0("log2fold_", tp, "_C1"),
          "log2fold"
        )

        dex_effect_col <- .first_present(dex_tbl, dex_effect_candidates)

        dex_gene_summary <- .summarize_gene_metric(
          tbl = dex_tbl,
          gene_col = "ensgene",
          sig_col = "padj",
          effect_col = dex_effect_col
        )
        dex_gene_summary <- dex_gene_summary[dex_gene_summary$ensgene %in% candidate_genes, , drop = FALSE]
        colnames(dex_gene_summary)[2:3] <- c("dex_best_padj", "dex_best_effect")
      })
      message("", tp, ": DEXSeq summary time = ", round(dex_time["elapsed"], 2), " sec")

      # DRIMSeq: use transcript-level table
      message("", tp, ": summarizing DRIMSeq transcript-level metrics")
      drim_time <- system.time({
        drim_tbl <- drimseq_results$results$tx_full_by_tp[[tp]]
        drim_effect_col <- .first_present(
          drim_tbl,
          c("proportion", "coef", "effect_size")
        )

        drim_gene_summary <- .summarize_gene_metric(
          tbl = drim_tbl,
          gene_col = "ensgene",
          sig_col = "padj",
          effect_col = drim_effect_col
        )
        drim_gene_summary <- drim_gene_summary[drim_gene_summary$ensgene %in% candidate_genes, , drop = FALSE]
        colnames(drim_gene_summary)[2:3] <- c("drim_best_padj", "drim_best_effect")
      })
      message("", tp, ": DRIMseq summary time = ", round(drim_time["elapsed"], 2), " sec")

      # SUPPA2: use event-level table
      message("", tp, ": summarizing SUPPA2 event-level metrics")
      suppa_time <- system.time({
        suppa_tbl <- suppa_results$results$event_full_by_tp[[tp]]

        suppa_gene_summary <- .summarize_gene_metric(
          tbl = suppa_tbl,
          gene_col = "ensgene",
          sig_col = "pvalue",
          effect_col = "dpsi"
        )
        suppa_gene_summary <- suppa_gene_summary[suppa_gene_summary$ensgene %in% candidate_genes, , drop = FALSE]
        colnames(suppa_gene_summary)[2:3] <- c("suppa_best_pvalue", "suppa_best_dpsi")
      })
      message("", tp, ": SUPPA2 summary time = ", round(suppa_time["elapsed"], 2), " sec")

      message("", tp, ": collapsing SUPPA2 event types by gene")
      suppa_event_types <- if (!is.null(suppa_tbl) && nrow(suppa_tbl) > 0) {
        tmp <- aggregate(
          event_type ~ ensgene,
          data = suppa_tbl[!is.na(suppa_tbl$event_type), c("ensgene", "event_type"), drop = FALSE],
          FUN = function(x) paste(sort(unique(as.character(x))), collapse = ",")
        )
        colnames(tmp)[2] <- "suppa_event_types"
        tmp
      } else {
        data.frame(ensgene = character(), suppa_event_types = character(), stringsAsFactors = FALSE)
      }

      message("", tp, ": annotating and merging candidate summaries")
      df <- data.frame(
        timepoint = tp,
        ensgene = candidate_genes,
        stringsAsFactors = FALSE
      )

      df <- merge(df, .annotate_ensgenes(df$ensgene, annot_df = annot_df), by = "ensgene", all.x = TRUE, sort = FALSE)
      df <- merge(df, dex_gene_summary, by = "ensgene", all.x = TRUE, sort = FALSE)
      df <- merge(df, drim_gene_summary, by = "ensgene", all.x = TRUE, sort = FALSE)
      df <- merge(df, suppa_gene_summary, by = "ensgene", all.x = TRUE, sort = FALSE)
      df <- merge(df, suppa_event_types, by = "ensgene", all.x = TRUE, sort = FALSE)


      message("", tp, ": scoring and ranking candidates")
      dex_sig_score <- 1 - .rescale(df$dex_best_padj)
      drim_sig_score <- 1 - .rescale(df$drim_best_padj)
      suppa_sig_score <- 1 - .rescale(df$suppa_best_pvalue)

      df$mean_sig_score <- rowMeans(
        cbind(dex_sig_score, drim_sig_score, suppa_sig_score),
        na.rm = TRUE
      )

      dex_eff_score <- .rescale(df$dex_best_effect)
      drim_eff_score <- .rescale(df$drim_best_effect)
      suppa_eff_score <- .rescale(df$suppa_best_dpsi)

      df$mean_effect_score <- rowMeans(
        cbind(dex_eff_score, drim_eff_score, suppa_eff_score),
        na.rm = TRUE
      )

      df$candidate_score <- rowMeans(
        cbind(df$mean_sig_score, df$mean_effect_score),
        na.rm = TRUE
      )

      df <- df[order(-df$candidate_score, -df$mean_effect_score, df$suppa_best_pvalue), , drop = FALSE]
      rownames(df) <- NULL
      message("", tp, ": done")
      df
    }),
    timepoints
  )

  candidate_summary_df <- do.call(rbind, lapply(timepoints, function(tp) {
    data.frame(
      timepoint = tp,
      n_deg = length(deg_genes_by_tp[[tp]]),
      n_dex = length(dex_genes_by_tp[[tp]]),
      n_drim = length(drim_genes_by_tp[[tp]]),
      n_suppa = length(suppa_genes_by_tp[[tp]]),
      n_common_all_3 = length(common_splicing_genes_by_tp[[tp]]),
      n_common_splicing_only = length(splicing_only_common_genes_by_tp[[tp]]),
      n_common_with_deg = length(common_with_deg_genes_by_tp[[tp]]),
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(focus_timepoint)) {
    focus_metric <- switch(focus_mode,
      splicing_only = "n_common_splicing_only",
      shared_with_deg = "n_common_with_deg"
    )

    focus_timepoint <- candidate_summary_df$timepoint[
      which.max(candidate_summary_df[[focus_metric]])
    ]
  }

  out <- list(
    meta = list(
      outdir = outdir,
      timepoints = timepoints,
      focus_mode = focus_mode,
      focus_timepoint = focus_timepoint
    ),
    tables = list(
      candidate_summary = candidate_summary_df,
      candidate_rankings_by_tp = candidate_rankings_by_tp,
      focus_candidates = candidate_rankings_by_tp[[focus_timepoint]]
    ),
    gene_sets = list(
      deg_genes_by_tp = deg_genes_by_tp,
      dex_genes_by_tp = dex_genes_by_tp,
      drim_genes_by_tp = drim_genes_by_tp,
      suppa_genes_by_tp = suppa_genes_by_tp,
      common_splicing_genes_by_tp = common_splicing_genes_by_tp,
      splicing_only_common_genes_by_tp = splicing_only_common_genes_by_tp,
      common_with_deg_genes_by_tp = common_with_deg_genes_by_tp,
      candidate_gene_sets_by_tp = candidate_gene_sets_by_tp
    )
  )

  rds_file <- file.path(
    outdir,
    paste0("fig_3_candidate_gene_summary_", focus_mode, ".rds")
  )

  saveRDS(out, file = rds_file)

  message("Saved analysis object: ", rds_file)
  message("Analysis complete")
  message("focus_timepoint = ", focus_timepoint)

  out
}

# --------------------------------------------------
# Supplemental Figure 2A - Candidate Gene Splicejam Plot
# --------------------------------------------------

plot_candidate_gene_splicejam <- function(
  gene,
  tp = "H3",
  suppa_results,
  conditions = c("C1", tp),
  event_id = NULL,
  zoom_region = NULL,
  pad = 250,
  minJunctionScore = 1,
  txdb_info = NULL,
  tx2geneDF = NULL,
  se = NULL,
  full_rel_heights = c(1, 2),
  all_rel_heights = c(3, 2),
  save_plot = TRUE
) {
  .check_pkg(c(
    "cowplot",
    "ggplot2"
  ))

  sashimi_data <- .prep_sashimi_data(
    gene = gene,
    conditions = conditions,
    txdb_info = txdb_info,
    tx2geneDF = tx2geneDF,
    se = se
  )

  zoom_info <- .get_zoom_region(
    suppa_results = suppa_results,
    gene = gene,
    tp = tp,
    event_id = event_id,
    pad = pad
  )

  if (is.null(zoom_region)) {
    zoom_region <- zoom_info$zoom_region
  }

  message("Selected SUPPA event for ", gene, " (", tp, "):")
  message("  event_id   = ", zoom_info$event_id)
  if (!is.na(zoom_info$event_type)) {
    message("  event_type = ", zoom_info$event_type)
  }
  message("  dPSI       = ", signif(zoom_info$dpsi, 4))
  message("  zoom_region = ", paste(zoom_region, collapse = " - "))

  p_full <- .plot_candidate_gene_splicejam_model_and_sashimi(
    sashimi_data = sashimi_data,
    tp = tp,
    conditions = conditions,
    zoom_region = zoom_region,
    minJunctionScore = minJunctionScore,
    show_x_title = FALSE,
    rel_heights = full_rel_heights,
    save_plot = FALSE
  )

  p_zoom <- .plot_candidate_gene_splicejam_sashimi_zoom_only(
    sashimi_data = sashimi_data,
    tp = tp,
    conditions = conditions,
    zoom_region = zoom_region,
    minJunctionScore = minJunctionScore,
    save_plot = FALSE
  )

  p_psi <- .plot_suppa_psi(
    suppa_results = suppa_results,
    gene = gene,
    tp = tp,
    conditions = conditions,
    event_id = zoom_info$event_id,
    zoom_region = zoom_region,
    show_title = FALSE,
    save_plot = FALSE
  )

  p_bottom <- cowplot::plot_grid(
    p_zoom$plot,
    p_psi$plot,
    nrow = 1,
    rel_widths = c(1, 1),
    align = "h",
    axis = "tb"
  )

  p_all <- cowplot::plot_grid(
    p_full$plot,
    p_bottom,
    ncol = 1,
    align = "v",
    axis = "lr",
    rel_heights = all_rel_heights
  )

  saved_file <- NULL
  if (isTRUE(save_plot)) {
    outdir <- .set_outdir(subdir = "fig_3")
    gene_dir <- file.path(outdir, gene)
    dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

    outfile <- file.path(gene_dir, paste0("fig_3_A_candidate_gene_splicejam_", gene, ".pdf"))
    ggplot2::ggsave(outfile, p_all, width = 12, height = 15)
    saved_file <- outfile
  }

  list(
    sashimi_data = sashimi_data,
    zoom_info = zoom_info,
    zoom_region = zoom_region,
    selected_event_id = zoom_info$event_id,
    plots = list(
      full = p_full,
      zoom = p_zoom,
      psi = p_psi,
      all = p_all
    ),
    saved_file = saved_file
  )
}

# --------------------------------------------------
# Supplemental Figure 2B - DRIMseq Ribbon Plot
# --------------------------------------------------

plot_candidate_gene_ribbon <- function(
  drimseq_results,
  gene,
  tp,
  conditions,
  fdr = 0.10,
  txdb_info = NULL,
  tx2geneDF = NULL,
  se = NULL,
  save_plot = TRUE,
  plot_width = 9,
  plot_height = 6
) {
  .check_pkg(c(
    "reshape2",
    "ggplot2",
    "splicejam",
    "SummarizedExperiment"
  ))

  .require_two_conditions(conditions, "plot_candidate_gene_ribbon")

  if (is.null(txdb_info)) {
    paths <- .get_runtime_paths()
    env_file <- paths$env_file
    resources_dir <- paths$resources_dir

    txdb_info <- .get_txdb(
      env_file = env_file,
      resources_dir = resources_dir,
      overwrite = FALSE
    )
  }

  if (is.null(tx2geneDF)) {
    tx2geneDF <- .make_tx2gene_df(txdb_info)
  }

  # Resolve gene_id from symbol
  gene_id <- unique(as.character(tx2geneDF$gene_id[tx2geneDF$gene_name == gene]))
  gene_id <- gene_id[nzchar(gene_id)]

  if (length(gene_id) == 0) {
    stop("Could not resolve gene symbol to Ensembl ID: ", gene)
  }
  if (length(gene_id) > 1) {
    stop("Gene symbol maps to multiple Ensembl IDs: ", gene)
  }

  gene_id_unv <- .strip_ens_version(gene_id)

  # Retrieve expressed transcripts
  det <- .get_detected_transcripts_for_gene(
    gene = gene_id,
    conditions = conditions,
    tx2geneDF = tx2geneDF,
    se = se,
    force = FALSE
  )

  detected_tx <- det$detectedTx
  gene_symbol <- det$gene_name

  if (length(detected_tx) == 0) {
    stop("No detected transcripts found for gene: ", gene_id)
  }

  dm_obj <- drimseq_results$by_tp[[tp]]$dm
  if (is.null(dm_obj)) {
    stop("Could not find DRIMSeq dm object for timepoint: ", tp)
  }

  drim_tx_tbl <- drimseq_results$results$tx_full_by_tp[[tp]]
  if (is.null(drim_tx_tbl) || !is.data.frame(drim_tx_tbl)) {
    stop("Could not find transcript-level DRIMSeq results table for timepoint: ", tp)
  }

  # Match gene from the DRIMSeq results table
  gene_rows <- drim_tx_tbl[
    .strip_ens_version(drim_tx_tbl$gene_id_full) == gene_id_unv, ,
    drop = FALSE
  ]

  if (nrow(gene_rows) == 0) {
    stop(
      "Could not find gene in DRIMSeq transcript results table.\n",
      "Input gene_id: ", gene_id, "\n",
      "Unversioned gene_id: ", gene_id_unv
    )
  }

  gene_id_drim <- unique(as.character(gene_rows$gene_id_full))[1]
  fit_full <- as.matrix(dm_obj@fit_full[[gene_id_drim]])

  if (nrow(fit_full) == 0) {
    stop(
      "Resolved DRIMSeq gene key but fit/count matrices are empty.\n",
      "Resolved gene_id_full: ", gene_id_drim
    )
  }

  # Transcript IDs come from the DRIMSeq results table
  tx_id_col <- "feature_id"
  if (!tx_id_col %in% colnames(gene_rows)) {
    stop("feature_id column not found in DRIMSeq transcript results.")
  }

  feature_ids <- unique(as.character(gene_rows[[tx_id_col]]))
  feature_ids <- feature_ids[nzchar(feature_ids)]

  if (length(feature_ids) != nrow(fit_full)) {
    stop(
      "Transcript count mismatch between DRIMSeq results table and fit matrix.\n",
      "length(feature_ids) = ", length(feature_ids), "\n",
      "nrow(fit_full) = ", nrow(fit_full)
    )
  }

  proportions_fit <- fit_full
  proportions_fit[!is.finite(proportions_fit)] <- NA_real_

  prop_fit <- as.data.frame(
    proportions_fit,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  prop_fit$feature_id <- feature_ids
  prop_fit <- prop_fit[, c("feature_id", setdiff(colnames(prop_fit), "feature_id")), drop = FALSE]

  prop_fit_long <- reshape2::melt(
    prop_fit,
    id.vars = "feature_id",
    variable.name = "sample_id",
    value.name = "proportion"
  )

  infer_condition <- function(x) sub("_.*$", "", x)
  prop_fit_long$group <- factor(infer_condition(prop_fit_long$sample_id), levels = conditions)

  # One fitted proportion per transcript per condition
  prop_fit_long <- aggregate(
    proportion ~ feature_id + group,
    data = prop_fit_long,
    FUN = median
  )

  # Keep expressed transcripts
  keep_tx <- intersect(detected_tx, unique(as.character(prop_fit_long$feature_id)))
  if (length(keep_tx) == 0) {
    stop("None of the detected transcripts are present in the DRIMSeq fit for ", gene_id)
  }

  # Ensure all transcripts are present in all groups
  all_combos <- expand.grid(
    feature_id = keep_tx,
    group = conditions,
    stringsAsFactors = FALSE
  )

  prop_fit_long <- merge(
    all_combos,
    prop_fit_long,
    by = c("feature_id", "group"),
    all.x = TRUE
  )

  # Fill missing proportions with 0
  prop_fit_long$proportion[is.na(prop_fit_long$proportion)] <- 0

  # Enforce consistent ordering
  prop_fit_long$feature_id <- factor(prop_fit_long$feature_id, levels = keep_tx)
  prop_fit_long$group <- factor(prop_fit_long$group, levels = conditions)

  # Label significant transcripts
  sig_tx <- gene_rows[
    !is.na(gene_rows$padj) & gene_rows$padj < fdr,
    tx_id_col,
    drop = TRUE
  ]
  sig_tx <- intersect(as.character(sig_tx), keep_tx)

  feature_palette <- unname(category_base_colors[seq_len(length(keep_tx))])
  names(feature_palette) <- keep_tx

  feature_labels <- stats::setNames(
    ifelse(keep_tx %in% sig_tx, paste0(keep_tx, " *"), keep_tx),
    keep_tx
  )

  condition_labels <- .get_condition_labels(conditions)

  p <- .plot_drimseq_ribbon(
    prop_fit = prop_fit_long,
    main = gene_symbol,
    feature_colors = feature_palette
  ) +
    ggplot2::scale_fill_manual(
      values = feature_palette,
      breaks = keep_tx,
      labels = feature_labels,
      name = "Transcripts"
    ) +
    ggplot2::scale_x_discrete(
      labels = condition_labels
    ) +
    ggplot2::labs(
      subtitle = "DRIMSeq Differential Transcript Usage"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 18),
      plot.subtitle = ggplot2::element_text(size = 12.5),
      axis.title.x = ggplot2::element_text(size = 15),
      axis.title.y = ggplot2::element_text(size = 15, margin = ggplot2::margin(r = 10)),
      axis.text.x = ggplot2::element_text(size = 15, margin = ggplot2::margin(t = 10)),
      legend.title = ggplot2::element_text(size = 15),
      legend.text = ggplot2::element_text(size = 12.5)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        keyheight = grid::unit(0.75, "cm"),
        keywidth = grid::unit(0.75, "cm")
      )
    )

  if (isTRUE(save_plot)) {
    outdir <- .set_outdir(subdir = "fig_3")
    gene_dir <- file.path(outdir, gene)
    dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

    outfile <- file.path(gene_dir, paste0("fig_3_B_candidate_gene_ribbon_", gene_symbol, ".pdf"))
    ggplot2::ggsave(outfile, p, width = plot_width, height = plot_height)
  }

  p
}

# --------------------------------------------------
# Supplemental Figure 2C - DEXseq DEU Plot
# --------------------------------------------------

plot_candidate_gene_deu <- function(
  dexseq_results,
  gene,
  tp,
  conditions = c("C1", tp),
  txdb_info = NULL,
  tx2geneDF = NULL,
  color = NULL,
  fdr = 0.10,
  fitExpToVar = "condition",
  norCounts = FALSE,
  expression = FALSE,
  splicing = TRUE,
  displayTranscripts = FALSE,
  names = FALSE,
  legend = FALSE,
  cex.axis = 1,
  cex = 1,
  lwd = 3,
  save_plot = TRUE,
  plot_width = 15,
  plot_height = 9
) {
  .check_pkg("DEXSeq")

  paths <- .get_runtime_paths()

  if (is.null(txdb_info)) {
    paths <- .get_runtime_paths()
    env_file <- paths$env_file
    resources_dir <- paths$resources_dir

    txdb_info <- .get_txdb(
      env_file = env_file,
      resources_dir = resources_dir,
      overwrite = FALSE
    )
  }

  if (is.null(tx2geneDF)) {
    tx2geneDF <- .make_tx2gene_df(txdb_info)
  }

  condition_labels <- .get_condition_labels(conditions)

  if (is.null(color)) {
    color <- unname(.get_condition_shaded_colors(conditions, shade = "base"))
  }

  # Resolve gene_id from symbol
  gene_id <- unique(as.character(tx2geneDF$gene_id[tx2geneDF$gene_name == gene]))
  gene_id <- gene_id[nzchar(gene_id)]

  if (length(gene_id) == 0) {
    stop("Could not resolve gene symbol to Ensembl ID: ", gene)
  }
  if (length(gene_id) > 1) {
    stop("Gene symbol maps to multiple Ensembl IDs: ", gene)
  }

  if (!tp %in% names(dexseq_results$by_tp)) {
    stop("Timepoint not found in dexseq_results$by_tp: ", tp)
  }

  dxr <- dexseq_results$by_tp[[tp]]$dxr
  if (is.null(dxr)) {
    stop("DEXSeq results object 'dxr' not found for timepoint: ", tp)
  }

  if (isTRUE(save_plot)) {
    outdir <- .set_outdir(subdir = "fig_3")
    gene_dir <- file.path(outdir, gene)
    dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

    outfile <- file.path(gene_dir, paste0("fig_3_C_candidate_gene_deu_", gene, ".pdf"))
    grDevices::pdf(outfile, width = plot_width, height = plot_height)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  DEXSeq::plotDEXSeq(
    dxr,
    gene_id,
    FDR = fdr,
    fitExpToVar = fitExpToVar,
    norCounts = norCounts,
    expression = expression,
    splicing = splicing,
    displayTranscripts = displayTranscripts,
    names = names,
    legend = legend,
    color = color,
    cex.axis = cex.axis,
    cex = cex,
    lwd = lwd
  )

  # Put custom legend in top panel
  graphics::par(mfg = c(1, 1))

  usr <- graphics::par("usr")

  graphics::legend(
    x = usr[2] - 0.20 * diff(usr[1:2]),
    y = usr[4] - 0.10 * diff(usr[3:4]),
    legend = condition_labels,
    col = color,
    lwd = 3,
    pt.cex = 1.5,
    cex = 1.5,
    bty = "n",
    xjust = 1,
    yjust = 1
  )

  plotted <- grDevices::recordPlot()

  if (isTRUE(save_plot)) {
    attr(plotted, "outfile") <- outfile
  }

  plotted
}

# --------------------------------------------------
# Convenience Wrapper
# --------------------------------------------------

plot_fig_3_all <- function(
  gene,
  tp,
  conditions = c("C1", tp),
  event_id = NULL,
  zoom_region = NULL,
  dexseq_results,
  drimseq_results,
  suppa_results,
  save_plot = TRUE
) {
  stopifnot(is.character(gene), length(gene) == 1, nzchar(gene))
  stopifnot(is.character(tp), length(tp) == 1, nzchar(tp))
  .require_two_conditions(conditions, "plot_fig_3_all")

  message("Generating Figure 3 Plots for ", gene)
  message("Timepoint: ", tp)
  message("Conditions: ", paste(conditions, collapse = ", "))

  paths <- .get_runtime_paths()

  txdb_info <- .get_txdb(
    env_file = paths$env_file,
    resources_dir = paths$resources_dir,
    overwrite = FALSE
  )

  tx2geneDF <- .make_tx2gene_df(txdb_info)

  se <- .get_se(
    coldata_tsv = paths$coldata_tsv,
    env_file = paths$env_file,
    resources_dir = paths$resources_dir,
    indexDir = paths$indexDir,
    force = FALSE
  )

  # Panel A - Splicejam + PSI
  message("[A] Building splicejam plot with PSI included")
  splicejam_call <- list(
    gene = gene,
    tp = tp,
    conditions = conditions,
    event_id = event_id,
    zoom_region = zoom_region,
    suppa_results = suppa_results,
    txdb_info = txdb_info,
    tx2geneDF = tx2geneDF,
    se = se,
    save_plot = save_plot
  )
  p_splicejam <- do.call(plot_candidate_gene_splicejam, splicejam_call)

  # Panel B - DRIMSeq Ribbon
  message("[B] Building DRIMSeq DTU plot")
  ribbon_call <- c(
    list(
      drimseq_results = drimseq_results,
      gene = gene,
      tp = tp,
      conditions = conditions,
      txdb_info = txdb_info,
      tx2geneDF = tx2geneDF,
      se = se,
      save_plot = save_plot
    )
  )
  p_ribbon <- do.call(plot_candidate_gene_ribbon, ribbon_call)

  # Panel C - DEXSeq DEU
  message("[C] Building DEXSeq DEU plot")
  deu_call <- c(
    list(
      dexseq_results = dexseq_results,
      gene = gene,
      tp = tp,
      conditions = conditions,
      txdb_info = txdb_info,
      tx2geneDF = tx2geneDF,
      save_plot = save_plot
    )
  )
  p_deu <- do.call(plot_candidate_gene_deu, deu_call)

  out <- list(
    gene = gene,
    tp = tp,
    conditions = conditions,
    splicejam = p_splicejam,
    ribbon = p_ribbon,
    deu = p_deu
  )

  message("Done: ", gene, " (", tp, ")")

  invisible(out)
}
