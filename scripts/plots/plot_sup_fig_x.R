# ============================================================
# plot_sup_fig_x_v2.R
#
# Supplemental analyses addressing the 1 h hypoxia comparison:
#   A. Gene-level sample PCA
#   B. Transcript-level sample PCA
#   C. HALLMARK_HYPOXIA GSEA
#   D. Exhaustive C1/H1 balanced-label permutation analysis
#   E. Exact balanced-label test of sample-level SUPPA2 PSI separation
#   F. Upper-triangular pairwise SUPPA2 PSI difference matrix
#
# The displayed gene PCA uses all genes in the DESeq2 tested universe. The
# displayed transcript PCA uses log2(Salmon TPM + 0.5) for expressed
# transcripts whose parent genes are in that universe. Separate sensitivity
# plots compare the top 500, 1,000, and 5,000 variable features with the full
# retained universes.
#
# PCA sample labels are repelled without connector segments. The plotting
# wrapper writes both the standard PCA panels and alternate versions with
# lightly filled enclosing ellipses. The envelopes are calculated in normalized
# PCA space so they can include conditions represented by only two samples;
# they are visual clustering guides, not inferential confidence regions.
#
# The permutation workflow reuses cached count matrices. It does not recount
# Salmon or BAM files, but it refits label-dependent model coefficients because
# fitted p-values cannot be made valid for a new labeling by changing colData
# alone. By default, DRIMSeq reuses the cached observed H1 gene-wise precision
# estimates to avoid repeating its expensive per-gene grid search. This is a
# conditional/descriptive sensitivity analysis; use
# drimseq_precision_mode = "refit" for a fully re-estimated DRIMSeq permutation.
# SUPPA2 commands are launched through PROJECT_ROOT/scripts/run by default, so
# psiPerEvent and diffSplice execute in the project's dedicated suppa2 conda
# environment even when suppa.py is not visible on the R process PATH.
#
# For three C1 and three H1 samples, choose(6, 3) / 2 = 10 unique,
# label-symmetric 3-vs-3 partitions. One is observed, leaving exactly nine
# alternative null partitions. The default is therefore exhaustive.
# The permutation figure displays only the four individual analysis tools in a
# 2 x 2 grid. Every observed and alternative assignment is labeled with its
# number of significant genes, and the shared assignment key is placed directly
# below the tool panels.
#
# In the combined figure, the gene- and transcript-level PCA panels are stacked
# in the upper-left column and share a condition legend immediately below that
# column. The two GSEA panels are stacked in the upper-right column. The
# permutation grid and its key occupy the middle section, and the two PSI panels
# occupy the bottom section. The PSI heatmap uses a vertical legend on its right.
# ============================================================

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT")
if (!nzchar(PROJECT_ROOT)) stop("PROJECT_ROOT env var not set.")

HELPERS_FILE <- file.path(PROJECT_ROOT, "scripts/utils/helpers.R")
if (!file.exists(HELPERS_FILE)) stop("helpers.R not found at: ", HELPERS_FILE)
source(HELPERS_FILE)

COLORS_FILE <- file.path(PROJECT_ROOT, "scripts/utils/color_palette.R")
if (!file.exists(COLORS_FILE)) stop("color_palette.R not found at: ", COLORS_FILE)
source(COLORS_FILE)

.SUP_FIG_X_ANALYSIS_VERSION <- "2026-07-29-v8"
.SUP_PERMUTATION_CACHE_VERSION <- c(
  DESeq2 = "v3_full_timecourse_gene_tximeta",
  DEXSeq = "v1",
  DRIMSeq = "v2",
  SUPPA2 = "v1"
)

# --------------------------------------------------
# General helpers
# --------------------------------------------------

.sup_check_pkg <- function(pkgs) {
  missing <- pkgs[
    !vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
  ]
  if (length(missing) > 0) {
    stop(
      "Missing required package(s): ", paste(missing, collapse = ", "),
      "\nInstall them in the project R environment and rerun."
    )
  }
  invisible(TRUE)
}

.sup_strip_version <- function(x) {
  sub("\\..*$", "", as.character(x))
}

.sup_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
  invisible(path)
}

.sup_read_result <- function(x, default_path, label, required = TRUE) {
  if (is.list(x) || methods::is(x, "SummarizedExperiment")) {
    return(x)
  }

  path <- x %||% default_path
  if (is.null(path) || !file.exists(path)) {
    if (isTRUE(required)) {
      stop(label, " cache not found: ", path %||% "<no path supplied>")
    }
    return(NULL)
  }

  message("Loading cached ", label, ": ", path)
  readRDS(path)
}

.sup_find_dexseq_cache <- function(results_dir = .get_results_dir()) {
  dex_dir <- file.path(results_dir, "analysis", "dexseq")
  preferred <- file.path(
    dex_dir,
    "dexseq_results.C1.summarizeOverlaps.multiOverlap.all.rds"
  )
  if (file.exists(preferred)) return(preferred)

  candidates <- list.files(
    dex_dir,
    pattern = "^dexseq_results\\.C1\\..*\\.rds$",
    full.names = TRUE
  )
  if (length(candidates) == 0) return(preferred)

  summarize_candidates <- candidates[
    grepl("summarizeOverlaps", basename(candidates), fixed = TRUE)
  ]
  if (length(summarize_candidates) > 0) {
    candidates <- summarize_candidates
  }

  candidates <- sort(candidates)
  if (length(candidates) > 1) {
    message(
      "Multiple DEXSeq caches found; using: ", candidates[1],
      "\nPass dexseq_results explicitly to choose another."
    )
  }
  candidates[1]
}

.sup_condition_labels <- c(
  C1 = "Normoxia",
  H1 = "Hypoxia (1 h)",
  H3 = "Hypoxia (3 h)",
  H24 = "Hypoxia (24 h)"
)

.sup_condition_colors <- function() {
  fallback <- c(
    C1 = "#222222",
    H1 = "#0072B2",
    H3 = "#D55E00",
    H24 = "#CC79A7"
  )

  if (!exists("timepoint_base_colors", inherits = TRUE)) {
    return(fallback)
  }

  palette_object <- get("timepoint_base_colors", inherits = TRUE)
  palette_keys <- c(
    C1 = "Normoxia",
    H1 = "Hypoxia (1H)",
    H3 = "Hypoxia (3H)",
    H24 = "Hypoxia (24H)"
  )
  if (!all(palette_keys %in% names(palette_object))) {
    return(fallback)
  }

  stats::setNames(unname(palette_object[palette_keys]), names(palette_keys))
}

.sup_tool_colors <- function() {
  fallback <- c(
    DESeq2 = "#4C78A8",
    DEXSeq = "#F58518",
    DRIMSeq = "#54A24B",
    SUPPA2 = "#E45756"
  )

  if (exists("splicing_tool_colors", inherits = TRUE)) {
    pal <- get("splicing_tool_colors", inherits = TRUE)
    for (nm in intersect(c("DEXSeq", "DRIMSeq", "SUPPA2"), names(pal))) {
      fallback[nm] <- unname(pal[nm])
    }
  }
  fallback
}

.sup_analysis_assay <- function(se, preferred) {
  available <- SummarizedExperiment::assayNames(se)
  selected <- preferred[preferred %in% available][1]
  if (is.na(selected)) {
    stop(
      "None of the requested assays were found. Requested: ",
      paste(preferred, collapse = ", "),
      "; available: ", paste(available, collapse = ", ")
    )
  }
  SummarizedExperiment::assay(se, selected)
}

.sup_validate_conditions <- function(coldata, required = c("C1", "H1")) {
  if (!("condition" %in% colnames(coldata))) {
    stop("Expected a 'condition' column in sample metadata.")
  }
  observed <- unique(as.character(coldata$condition))
  missing <- setdiff(required, observed)
  if (length(missing) > 0) {
    stop("Missing required condition(s): ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

# --------------------------------------------------
# PCA analysis
# --------------------------------------------------

.sup_get_tested_gene_universe <- function(deseq_results) {
  tested <- deseq_results$results$tested_universe_all
  if (is.null(tested)) {
    tested_by_tp <- deseq_results$results$tested_universe_by_tp
    if (is.null(tested_by_tp)) {
      stop(
        "DESeq2 tested universe not found. Expected ",
        "deseq_results$results$tested_universe_all or tested_universe_by_tp."
      )
    }
    tested <- unique(unlist(tested_by_tp, use.names = FALSE))
  }

  tested <- unique(.sup_strip_version(tested))
  tested[!is.na(tested) & nzchar(tested)]
}

.sup_resolve_tx2gene_columns <- function(tx2gene) {
  tx_map <- as.data.frame(tx2gene, stringsAsFactors = FALSE)
  tx_col <- c("TXNAME", "transcript_id", "tx_id")
  gene_col <- c("GENEID", "gene_id", "ensgene")
  tx_col <- tx_col[tx_col %in% colnames(tx_map)][1]
  gene_col <- gene_col[gene_col %in% colnames(tx_map)][1]
  if (is.na(tx_col) || is.na(gene_col)) {
    stop("tx2gene must contain transcript and gene ID columns.")
  }

  list(data = tx_map, tx_col = tx_col, gene_col = gene_col)
}

.sup_prcomp_scores <- function(
  mat,
  coldata,
  ntop = NULL,
  feature_set,
  selection_label
) {
  mat <- as.matrix(mat)
  mat[!is.finite(mat)] <- NA_real_
  keep <- rowSums(is.finite(mat)) == ncol(mat)
  mat <- mat[keep, , drop = FALSE]

  row_var <- apply(mat, 1, stats::var)
  keep <- is.finite(row_var) & row_var > 0
  mat <- mat[keep, , drop = FALSE]
  row_var <- row_var[keep]

  if (nrow(mat) < 2) stop("Too few variable features for PCA.")

  if (is.null(ntop)) {
    selected <- seq_len(nrow(mat))
  } else {
    ntop <- min(as.integer(ntop), nrow(mat))
    selected <- order(row_var, decreasing = TRUE)[seq_len(ntop)]
  }
  pca <- stats::prcomp(
    t(mat[selected, , drop = FALSE]),
    center = TRUE,
    scale. = FALSE
  )

  scores <- as.data.frame(pca$x[, seq_len(min(4L, ncol(pca$x))), drop = FALSE])
  scores$sample_id <- rownames(scores)
  idx <- match(scores$sample_id, rownames(coldata))
  if (anyNA(idx)) stop("PCA sample names do not match sample metadata.")
  scores$condition <- as.character(coldata$condition[idx])
  scores$condition_label <- unname(.sup_condition_labels[scores$condition])
  scores$condition_label[is.na(scores$condition_label)] <- scores$condition[
    is.na(scores$condition_label)
  ]

  variance_explained <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

  list(
    pca = pca,
    scores = scores,
    variance_explained = variance_explained,
    n_features = length(selected),
    n_available_features = nrow(mat),
    feature_set = feature_set,
    selection_label = selection_label,
    selected_features = rownames(mat)[selected]
  )
}

.sup_run_pca_grid <- function(
  mat,
  coldata,
  feature_numbers,
  all_key,
  all_label,
  feature_noun
) {
  feature_numbers <- sort(unique(as.integer(feature_numbers)))
  feature_numbers <- feature_numbers[
    is.finite(feature_numbers) & feature_numbers > 0
  ]

  top_results <- stats::setNames(
    lapply(feature_numbers, function(n) {
      .sup_prcomp_scores(
        mat = mat,
        coldata = coldata,
        ntop = n,
        feature_set = paste0("top_", n),
        selection_label = paste0(
          "Top ", format(n, big.mark = ","), " variable ", feature_noun
        )
      )
    }),
    paste0("top_", feature_numbers)
  )

  all_result <- .sup_prcomp_scores(
    mat = mat,
    coldata = coldata,
    ntop = NULL,
    feature_set = all_key,
    selection_label = all_label
  )

  c(top_results, stats::setNames(list(all_result), all_key))
}

.sup_run_gene_pca <- function(
  deseq_results,
  feature_numbers = c(500L, 1000L, 5000L)
) {
  .sup_check_pkg(c("DESeq2", "SummarizedExperiment"))
  dds <- deseq_results$dds
  if (is.null(dds)) stop("Expected deseq_results$dds.")

  coldata <- as.data.frame(SummarizedExperiment::colData(dds))
  .sup_validate_conditions(coldata)

  vsd <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
  mat <- SummarizedExperiment::assay(vsd)
  tested_genes <- .sup_get_tested_gene_universe(deseq_results)
  keep <- .sup_strip_version(rownames(mat)) %in% tested_genes
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) < 2) {
    stop("Fewer than two DESeq2-tested genes were found in the VST matrix.")
  }

  list(
    results = .sup_run_pca_grid(
      mat = mat,
      coldata = coldata,
      feature_numbers = feature_numbers,
      all_key = "all_tested",
      all_label = "All DESeq2-tested genes",
      feature_noun = "genes"
    ),
    tested_gene_ids = tested_genes,
    transformation = "Blind variance-stabilizing transformation",
    universe_definition = paste(
      "Union of genes with non-missing DESeq2 adjusted p-values",
      "across H1, H3, and H24 contrasts."
    )
  )
}

.sup_run_transcript_pca <- function(
  transcript_se,
  tx2gene,
  tested_gene_ids,
  feature_numbers = c(500L, 1000L, 5000L),
  min_tpm = 1,
  min_samples = 3L,
  pseudocount = 0.5
) {
  .sup_check_pkg("SummarizedExperiment")

  if (!("abundance" %in% SummarizedExperiment::assayNames(transcript_se))) {
    stop(
      "Transcript PCA requires the Salmon/tximeta 'abundance' assay (TPM). ",
      "Raw estimated transcript counts are intentionally not used."
    )
  }
  abundance <- SummarizedExperiment::assay(transcript_se, "abundance")
  coldata <- as.data.frame(SummarizedExperiment::colData(transcript_se))
  .sup_validate_conditions(coldata)

  tx_info <- .sup_resolve_tx2gene_columns(tx2gene)
  tx_ids <- .sup_strip_version(rownames(abundance))
  map_idx <- match(
    tx_ids,
    .sup_strip_version(tx_info$data[[tx_info$tx_col]])
  )
  parent_gene <- .sup_strip_version(
    tx_info$data[[tx_info$gene_col]][map_idx]
  )

  expressed <- rowSums(abundance >= min_tpm, na.rm = TRUE) >= min_samples
  parent_tested <- !is.na(parent_gene) & parent_gene %in% tested_gene_ids
  keep <- expressed & parent_tested
  if (sum(keep) < 2) {
    stop(
      "Fewer than two transcripts passed the expression and tested-gene filters."
    )
  }

  mat <- log2(abundance[keep, , drop = FALSE] + pseudocount)

  list(
    results = .sup_run_pca_grid(
      mat = mat,
      coldata = coldata,
      feature_numbers = feature_numbers,
      all_key = "all_expressed",
      all_label = paste0(
        "All expressed transcripts (TPM >= ",
        format(min_tpm),
        " in >= ", min_samples, " samples)"
      ),
      feature_noun = "transcripts"
    ),
    expressed_transcript_ids = rownames(mat),
    parent_gene_ids = parent_gene[keep],
    transformation = paste0("log2(Salmon TPM + ", pseudocount, ")"),
    universe_definition = paste0(
      "Transcripts with TPM >= ", min_tpm, " in >= ", min_samples,
      " samples and a parent gene in the DESeq2-tested universe."
    ),
    min_tpm = min_tpm,
    min_samples = min_samples,
    pseudocount = pseudocount
  )
}

# --------------------------------------------------
# GSEA analysis
# --------------------------------------------------

.sup_reduce_rank <- function(values, ids) {
  ids <- .sup_strip_version(ids)
  keep <- is.finite(values) & !is.na(ids) & nzchar(ids)
  values <- as.numeric(values[keep])
  ids <- ids[keep]

  ord_abs <- order(abs(values), decreasing = TRUE)
  values <- values[ord_abs]
  ids <- ids[ord_abs]
  keep_unique <- !duplicated(ids)

  out <- values[keep_unique]
  names(out) <- ids[keep_unique]
  sort(out, decreasing = TRUE)
}

.sup_gene_rank <- function(deseq_results) {
  .sup_check_pkg(c("DESeq2", "SummarizedExperiment"))
  dds <- deseq_results$dds
  coldata <- as.data.frame(SummarizedExperiment::colData(dds))
  .sup_validate_conditions(coldata)

  res <- DESeq2::results(
    dds,
    contrast = c("condition", "H1", "C1"),
    independentFiltering = FALSE
  )
  .sup_reduce_rank(res$stat, rownames(res))
}

.sup_transcript_derived_gene_rank <- function(transcript_se, tx2gene) {
  .sup_check_pkg(c("limma", "SummarizedExperiment"))

  abundance <- .sup_analysis_assay(
    transcript_se,
    preferred = c("abundance", "counts")
  )
  coldata <- as.data.frame(SummarizedExperiment::colData(transcript_se))
  .sup_validate_conditions(coldata)

  keep_samples <- as.character(coldata$condition) %in% c("C1", "H1")
  sample_ids <- rownames(coldata)[keep_samples]
  abundance <- abundance[, sample_ids, drop = FALSE]
  condition <- factor(
    as.character(coldata$condition[keep_samples]),
    levels = c("C1", "H1")
  )

  if ("abundance" %in% SummarizedExperiment::assayNames(transcript_se)) {
    keep_tx <- rowSums(abundance >= 1, na.rm = TRUE) >= 3
    expr <- log2(abundance[keep_tx, , drop = FALSE] + 0.5)
  } else {
    keep_tx <- rowSums(abundance >= 10, na.rm = TRUE) >= 3
    expr <- log2(abundance[keep_tx, , drop = FALSE] + 1)
  }

  design <- stats::model.matrix(~0 + condition)
  colnames(design) <- levels(condition)
  fit <- limma::lmFit(expr, design)
  contrast <- limma::makeContrasts(H1 - C1, levels = design)
  fit <- limma::contrasts.fit(fit, contrast)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  tx_stat <- fit$t[, 1]

  tx_info <- .sup_resolve_tx2gene_columns(tx2gene)

  tx_ids <- .sup_strip_version(rownames(expr))
  map_idx <- match(
    tx_ids,
    .sup_strip_version(tx_info$data[[tx_info$tx_col]])
  )
  gene_ids <- tx_info$data[[tx_info$gene_col]][map_idx]
  .sup_reduce_rank(tx_stat, gene_ids)
}

.sup_get_hypoxia_genes <- function(rank_ids, hypoxia_genes = NULL) {
  rank_ids <- .sup_strip_version(rank_ids)

  if (!is.null(hypoxia_genes)) {
    ids <- intersect(rank_ids, .sup_strip_version(hypoxia_genes))
    if (length(ids) < 10) {
      stop(
        "The supplied hypoxia_genes overlap fewer than 10 ranked genes. ",
        "Supply Ensembl gene IDs matching the analysis annotation."
      )
    }
    return(list(
      ids = unique(ids),
      source = "user-supplied",
      name = "HALLMARK_HYPOXIA"
    ))
  }

  .sup_check_pkg("msigdbr")
  msig_formals <- names(formals(msigdbr::msigdbr))
  msig <- if ("collection" %in% msig_formals) {
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  } else {
    msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  }

  if (!("gs_name" %in% colnames(msig))) {
    stop("msigdbr output does not contain gs_name.")
  }
  msig <- msig[msig$gs_name == "HALLMARK_HYPOXIA", , drop = FALSE]
  if (nrow(msig) == 0) stop("HALLMARK_HYPOXIA was not found in msigdbr.")

  id_candidates <- c(
    "ensembl_gene",
    "db_ensembl_gene",
    "gene_symbol",
    "db_gene_symbol"
  )
  id_candidates <- id_candidates[id_candidates %in% colnames(msig)]
  overlaps <- lapply(id_candidates, function(nm) {
    intersect(rank_ids, .sup_strip_version(msig[[nm]]))
  })
  overlap_n <- vapply(overlaps, length, integer(1))
  if (length(overlap_n) == 0 || max(overlap_n) < 10) {
    stop(
      "Could not match at least 10 HALLMARK_HYPOXIA genes to the ranked ",
      "gene IDs. Pass hypoxia_genes as Ensembl IDs to override."
    )
  }

  best <- which.max(overlap_n)
  list(
    ids = unique(overlaps[[best]]),
    source = paste0("msigdbr:", id_candidates[best]),
    name = "HALLMARK_HYPOXIA"
  )
}

.sup_running_enrichment <- function(ranks, genes, ranking_label) {
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  hits <- names(ranks) %in% genes
  n_hits <- sum(hits)
  n_miss <- length(hits) - n_hits
  if (n_hits == 0 || n_miss == 0) {
    stop("Cannot calculate an enrichment curve with zero hits or misses.")
  }

  hit_weights <- abs(ranks)
  hit_weights[!hits] <- 0
  running <- cumsum(hit_weights / sum(hit_weights)) -
    cumsum((!hits) / n_miss)

  data.frame(
    rank = seq_along(ranks),
    running_score = running,
    hit = hits,
    ranked_statistic = unname(ranks),
    gene_id = names(ranks),
    ranking = ranking_label,
    stringsAsFactors = FALSE
  )
}

.sup_run_one_gsea <- function(ranks, gene_set, ranking_label) {
  .sup_check_pkg("fgsea")
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  genes <- intersect(names(ranks), gene_set$ids)

  fg <- fgsea::fgseaMultilevel(
    pathways = stats::setNames(list(genes), gene_set$name),
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    eps = 0
  )
  fg <- as.data.frame(fg, stringsAsFactors = FALSE)
  fg$ranking <- ranking_label
  fg$gene_set_source <- gene_set$source

  list(
    result = fg,
    running = .sup_running_enrichment(ranks, genes, ranking_label),
    ranks = ranks,
    genes = genes
  )
}

.sup_run_hypoxia_gsea <- function(
  deseq_results,
  transcript_se,
  tx2gene,
  hypoxia_genes = NULL,
  include_transcript_ranking = TRUE
) {
  gene_rank <- .sup_gene_rank(deseq_results)
  gene_set <- .sup_get_hypoxia_genes(names(gene_rank), hypoxia_genes)

  out <- list(
    gene = .sup_run_one_gsea(
      ranks = gene_rank,
      gene_set = gene_set,
      ranking_label = "Gene-Level DESeq2 Ranking"
    )
  )

  if (isTRUE(include_transcript_ranking)) {
    tx_rank <- .sup_transcript_derived_gene_rank(transcript_se, tx2gene)
    tx_gene_set <- gene_set
    tx_gene_set$ids <- intersect(gene_set$ids, names(tx_rank))

    out$transcript_derived <- .sup_run_one_gsea(
      ranks = tx_rank,
      gene_set = tx_gene_set,
      ranking_label = "Transcript-Abundance-Derived Ranking"
    )
  }

  out$summary <- do.call(
    rbind,
    lapply(out[intersect(c("gene", "transcript_derived"), names(out))], `[[`, "result")
  )
  rownames(out$summary) <- NULL
  out$gene_set <- gene_set
  out
}

# --------------------------------------------------
# Pairwise sample-level SUPPA2 PSI analysis
# --------------------------------------------------

.sup_read_suppa_psi <- function(path) {
  if (!file.exists(path)) stop("SUPPA2 PSI file not found: ", path)

  psi <- utils::read.delim(
    path,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NA", "NaN", "nan", "-1")
  )
  psi <- as.matrix(psi)
  storage.mode(psi) <- "numeric"

  if (nrow(psi) == 0 || ncol(psi) == 0) {
    stop("SUPPA2 PSI file contains no event-by-sample values: ", path)
  }
  if (anyDuplicated(rownames(psi))) {
    stop("SUPPA2 PSI event identifiers are not unique in: ", path)
  }
  psi
}

.sup_run_pairwise_psi <- function(
  deseq_results,
  psi_dir = NULL,
  cutoff = 0.10
) {
  .sup_check_pkg("SummarizedExperiment")
  if (!is.numeric(cutoff) || length(cutoff) != 1L ||
      !is.finite(cutoff) || cutoff <= 0 || cutoff > 1) {
    stop("pairwise_psi_cutoff must be one finite number in (0, 1].")
  }

  psi_dir <- psi_dir %||% file.path(
    .get_results_dir(),
    "analysis", "suppa", "standard", "C1", "events"
  )
  c1_file <- file.path(psi_dir, "events_C1.psi")
  h1_file <- file.path(psi_dir, "events_H1.psi")
  c1_psi <- .sup_read_suppa_psi(c1_file)
  h1_psi <- .sup_read_suppa_psi(h1_file)

  common_events <- intersect(rownames(c1_psi), rownames(h1_psi))
  if (length(common_events) == 0) {
    stop("No common event identifiers were found in C1 and H1 PSI files.")
  }
  psi <- cbind(
    c1_psi[common_events, , drop = FALSE],
    h1_psi[common_events, , drop = FALSE]
  )
  if (anyDuplicated(colnames(psi))) {
    stop("Duplicate sample names were found across the C1 and H1 PSI files.")
  }

  coldata <- as.data.frame(
    SummarizedExperiment::colData(deseq_results$dds)
  )
  .sup_validate_conditions(coldata)
  keep_samples <- rownames(coldata)[
    as.character(coldata$condition) %in% c("C1", "H1")
  ]
  missing_samples <- setdiff(keep_samples, colnames(psi))
  if (length(missing_samples) > 0) {
    stop(
      "SUPPA2 PSI matrices are missing C1/H1 sample(s): ",
      paste(missing_samples, collapse = ", ")
    )
  }
  sample_condition <- setNames(
    as.character(coldata[keep_samples, "condition"]),
    keep_samples
  )
  sample_order <- c(
    keep_samples[sample_condition[keep_samples] == "C1"],
    keep_samples[sample_condition[keep_samples] == "H1"]
  )
  psi <- psi[, sample_order, drop = FALSE]
  complete_events <- rowSums(is.finite(psi)) == ncol(psi)
  n_complete_events <- sum(complete_events)
  if (n_complete_events == 0) {
    stop("No SUPPA2 events have finite PSI in all six C1/H1 samples.")
  }
  psi <- psi[complete_events, , drop = FALSE]

  pair_index <- utils::combn(sample_order, 2)
  pair_summary <- do.call(rbind, lapply(seq_len(ncol(pair_index)), function(i) {
    sample_1 <- pair_index[1, i]
    sample_2 <- pair_index[2, i]
    condition_1 <- unname(sample_condition[sample_1])
    condition_2 <- unname(sample_condition[sample_2])

    # Keep normoxia first in between-condition labels and tables.
    if (condition_1 == "H1" && condition_2 == "C1") {
      tmp <- sample_1
      sample_1 <- sample_2
      sample_2 <- tmp
      condition_1 <- "C1"
      condition_2 <- "H1"
    }

    delta <- psi[, sample_2] - psi[, sample_1]
    values <- abs(delta)
    pair_type <- if (condition_1 == "C1" && condition_2 == "C1") {
      "Normoxia-Normoxia"
    } else if (condition_1 == "H1" && condition_2 == "H1") {
      "Hypoxia (1 h)-Hypoxia (1 h)"
    } else {
      "Normoxia-Hypoxia (1 h)"
    }

    data.frame(
      sample_1 = sample_1,
      condition_1 = condition_1,
      sample_2 = sample_2,
      condition_2 = condition_2,
      pair_label = paste(sample_1, sample_2, sep = " vs "),
      pair_type = pair_type,
      n_events_tested = n_complete_events,
      n_events_abs_dpsi_ge_cutoff = sum(values >= cutoff),
      fraction_events_abs_dpsi_ge_cutoff = if (length(values) > 0) {
        mean(values >= cutoff)
      } else {
        NA_real_
      },
      median_abs_dpsi = if (length(values) > 0) {
        stats::median(values)
      } else {
        NA_real_
      },
      p95_abs_dpsi = if (length(values) > 0) {
        unname(stats::quantile(values, 0.95))
      } else {
        NA_real_
      },
      mean_abs_dpsi = if (length(values) > 0) mean(values) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  rownames(pair_summary) <- NULL

  heatmap <- expand.grid(
    sample_1 = sample_order,
    sample_2 = sample_order,
    stringsAsFactors = FALSE
  )
  heatmap$fraction_events_abs_dpsi_ge_cutoff <- vapply(
    seq_len(nrow(heatmap)),
    function(i) {
      a <- heatmap$sample_1[i]
      b <- heatmap$sample_2[i]
      if (a == b) return(0)
      values <- abs(psi[, b] - psi[, a])
      mean(values >= cutoff)
    },
    numeric(1)
  )
  heatmap$percent_events_abs_dpsi_ge_cutoff <-
    100 * heatmap$fraction_events_abs_dpsi_ge_cutoff
  heatmap$condition_1 <- unname(sample_condition[heatmap$sample_1])
  heatmap$condition_2 <- unname(sample_condition[heatmap$sample_2])
  heatmap$n_events_tested <- n_complete_events

  distance_matrix <- matrix(
    heatmap$fraction_events_abs_dpsi_ge_cutoff,
    nrow = length(sample_order),
    ncol = length(sample_order),
    dimnames = list(sample_order, sample_order)
  )

  split_info <- .sup_make_balanced_splits(deseq_results)
  all_splits <- c(list(observed = split_info$observed), split_info$null)
  membership <- .sup_membership_table(all_splits, split_info$coldata)
  split_descriptors <- .sup_split_descriptors(membership)
  partition_summary <- do.call(rbind, lapply(names(all_splits), function(split_id) {
    z <- all_splits[[split_id]]
    within_a <- utils::combn(match(z$A, sample_order), 2)
    within_b <- utils::combn(match(z$B, sample_order), 2)
    within_values <- c(
      distance_matrix[cbind(within_a[1, ], within_a[2, ])],
      distance_matrix[cbind(within_b[1, ], within_b[2, ])]
    )
    between_values <- as.vector(distance_matrix[z$A, z$B, drop = FALSE])
    descriptor <- split_descriptors[
      split_descriptors$split_id == split_id,
      ,
      drop = FALSE
    ]

    data.frame(
      split_id = split_id,
      split_short = descriptor$split_short,
      split_type = if (split_id == "observed") "Observed" else "Permuted",
      swap_label = descriptor$swap_label,
      mean_within_distance_percent = 100 * mean(within_values),
      mean_between_distance_percent = 100 * mean(between_values),
      between_minus_within_percent_points = 100 * (
        mean(between_values) - mean(within_values)
      ),
      stringsAsFactors = FALSE
    )
  }))
  rownames(partition_summary) <- NULL

  observed_statistic <- partition_summary$between_minus_within_percent_points[
    partition_summary$split_type == "Observed"
  ]
  if (length(observed_statistic) != 1L) {
    stop("Could not uniquely identify the observed PSI-distance partition.")
  }
  partition_summary$separation_rank <- rank(
    -partition_summary$between_minus_within_percent_points,
    ties.method = "min"
  )
  exact_p <- mean(
    partition_summary$between_minus_within_percent_points >=
      observed_statistic
  )
  partition_test <- data.frame(
    distance_definition = paste0(
      "Fraction of complete-case events with absolute delta PSI >= ",
      format(cutoff, nsmall = 2)
    ),
    n_events = n_complete_events,
    n_unique_balanced_partitions = nrow(partition_summary),
    observed_mean_within_distance_percent =
      partition_summary$mean_within_distance_percent[
        partition_summary$split_type == "Observed"
      ],
    observed_mean_between_distance_percent =
      partition_summary$mean_between_distance_percent[
        partition_summary$split_type == "Observed"
      ],
    observed_between_minus_within_percent_points = observed_statistic,
    observed_separation_rank =
      partition_summary$separation_rank[
        partition_summary$split_type == "Observed"
      ],
    exact_upper_tail_p = exact_p,
    stringsAsFactors = FALSE
  )

  list(
    pair_summary = pair_summary,
    heatmap = heatmap,
    partition_summary = partition_summary,
    partition_test = partition_test,
    sample_order = sample_order,
    sample_condition = sample_condition[sample_order],
    cutoff = cutoff,
    n_common_event_ids = length(common_events),
    n_complete_events = n_complete_events,
    n_common_events = n_complete_events,
    event_universe_definition = paste(
      "SUPPA2 event IDs shared by C1 and H1 PSI files with finite PSI",
      "in all six samples."
    ),
    psi_files = c(C1 = c1_file, H1 = h1_file)
  )
}

# --------------------------------------------------
# Balanced-label permutations
# --------------------------------------------------

.sup_make_balanced_splits <- function(deseq_results) {
  .sup_check_pkg("SummarizedExperiment")
  coldata <- as.data.frame(
    SummarizedExperiment::colData(deseq_results$dds)
  )
  .sup_validate_conditions(coldata)

  keep <- as.character(coldata$condition) %in% c("C1", "H1")
  c1h1 <- coldata[keep, , drop = FALSE]
  if (nrow(c1h1) != 6) {
    stop(
      "The exhaustive balanced design requires exactly six C1/H1 samples; ",
      "found ", nrow(c1h1), "."
    )
  }
  if (sum(as.character(c1h1$condition) == "C1") != 3 ||
      sum(as.character(c1h1$condition) == "H1") != 3) {
    stop("Expected three C1 and three H1 samples.")
  }

  sample_ids <- sort(rownames(c1h1))
  candidates <- utils::combn(sample_ids, 3, simplify = FALSE)
  canonical <- lapply(candidates, function(a) {
    b <- setdiff(sample_ids, a)
    key_a <- paste(sort(a), collapse = ",")
    key_b <- paste(sort(b), collapse = ",")
    if (key_a <= key_b) {
      list(A = sort(a), B = sort(b), key = paste(key_a, key_b, sep = "|"))
    } else {
      list(A = sort(b), B = sort(a), key = paste(key_b, key_a, sep = "|"))
    }
  })
  keys <- vapply(canonical, `[[`, character(1), "key")
  canonical <- canonical[!duplicated(keys)]
  canonical <- canonical[order(vapply(canonical, `[[`, character(1), "key"))]

  observed_c1 <- sort(rownames(c1h1)[as.character(c1h1$condition) == "C1"])
  observed_h1 <- sort(rownames(c1h1)[as.character(c1h1$condition) == "H1"])
  is_observed <- vapply(canonical, function(z) {
    setequal(z$A, observed_c1) || setequal(z$B, observed_c1)
  }, logical(1))
  if (sum(is_observed) != 1) stop("Could not uniquely identify the observed split.")

  observed <- list(
    A = observed_c1,
    B = observed_h1,
    key = paste(
      paste(observed_c1, collapse = ","),
      paste(observed_h1, collapse = ","),
      sep = "|"
    )
  )
  null <- canonical[!is_observed]
  names(null) <- sprintf("null_%02d", seq_along(null))

  list(
    observed = observed,
    null = null,
    coldata = c1h1,
    sample_ids = sample_ids
  )
}

.sup_select_null_splits <- function(null_splits, n_null_splits, seed) {
  n_available <- length(null_splits)
  if (is.infinite(n_null_splits) || n_null_splits >= n_available) {
    return(null_splits)
  }

  n_null_splits <- as.integer(n_null_splits)
  if (!is.finite(n_null_splits) || n_null_splits < 1) {
    stop("n_null_splits must be a positive integer or Inf.")
  }

  warning(
    "Running only ", n_null_splits, " of ", n_available,
    " alternative splits. This is a pilot/descriptive run, not an exhaustive ",
    "randomization analysis. A single shuffle should not be used as inferential evidence."
  )
  set.seed(seed)
  null_splits[sample(seq_along(null_splits), n_null_splits)]
}

.sup_membership_table <- function(splits, original_coldata) {
  do.call(rbind, lapply(names(splits), function(split_id) {
    z <- splits[[split_id]]
    sample_ids <- c(z$A, z$B)
    idx <- match(sample_ids, rownames(original_coldata))
    data.frame(
      split_id = split_id,
      sample_id = sample_ids,
      pseudo_group = rep(c("A", "B"), c(length(z$A), length(z$B))),
      original_condition = as.character(original_coldata$condition[idx]),
      stringsAsFactors = FALSE
    )
  }))
}

.sup_split_descriptors <- function(membership) {
  split_ids <- unique(membership$split_id)
  out <- do.call(rbind, lapply(split_ids, function(split_id) {
    z <- membership[membership$split_id == split_id, , drop = FALSE]
    n_c1_by_group <- tapply(
      z$original_condition == "C1",
      z$pseudo_group,
      sum
    )
    pseudo_c1_group <- names(n_c1_by_group)[which.max(n_c1_by_group)]
    pseudo_h1_group <- setdiff(unique(z$pseudo_group), pseudo_c1_group)

    pseudo_c1_members <- sort(z$sample_id[z$pseudo_group == pseudo_c1_group])
    pseudo_h1_members <- sort(z$sample_id[z$pseudo_group == pseudo_h1_group])
    swapped_c1 <- sort(z$sample_id[
      z$original_condition == "C1" &
        z$pseudo_group == pseudo_h1_group
    ])
    swapped_h1 <- sort(z$sample_id[
      z$original_condition == "H1" &
        z$pseudo_group == pseudo_c1_group
    ])
    split_short <- if (split_id == "observed") {
      "Observed"
    } else {
      paste0("P", as.integer(sub("^null_", "", split_id)))
    }
    swap_label <- if (split_id == "observed") {
      "Observed labels"
    } else {
      paste(
        paste(swapped_c1, collapse = ", "),
        paste(swapped_h1, collapse = ", "),
        sep = " <-> "
      )
    }

    data.frame(
      split_id = split_id,
      split_short = split_short,
      swap_label = swap_label,
      swapped_C1_sample = paste(swapped_c1, collapse = ", "),
      swapped_H1_sample = paste(swapped_h1, collapse = ", "),
      pseudo_C1_members = paste(pseudo_c1_members, collapse = ", "),
      pseudo_H1_members = paste(pseudo_h1_members, collapse = ", "),
      n_original_labels_preserved = sum(
        (z$original_condition == "C1" & z$pseudo_group == pseudo_c1_group) |
          (z$original_condition == "H1" & z$pseudo_group == pseudo_h1_group)
      ),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

.sup_summary_row <- function(
  split_id,
  tool,
  level,
  n_tested,
  n_significant,
  threshold,
  status = "ok"
) {
  data.frame(
    split_id = split_id,
    tool = tool,
    level = level,
    n_tested = as.integer(n_tested),
    n_significant = as.integer(n_significant),
    fraction_significant = if (
      is.finite(n_tested) && n_tested > 0
    ) {
      n_significant / n_tested
    } else {
      NA_real_
    },
    threshold = threshold,
    status = status,
    stringsAsFactors = FALSE
  )
}

.sup_run_deseq_split <- function(
  split_id,
  split,
  deseq_results,
  gene_se,
  alpha
) {
  .sup_check_pkg(c("DESeq2", "SummarizedExperiment", "S4Vectors"))
  source_dds <- deseq_results$dds
  input_coldata <- as.data.frame(SummarizedExperiment::colData(gene_se))

  if (!setequal(colnames(gene_se), colnames(source_dds))) {
    stop(
      "The cached gene-level tximeta object and manuscript DESeq2 object ",
      "do not contain the same samples."
    )
  }
  if (!all(c(split$A, split$B) %in% rownames(input_coldata))) {
    stop("Gene-level tximeta object is missing permutation samples.")
  }
  if (!("condition" %in% colnames(input_coldata))) {
    stop("Expected a 'condition' column in the gene-level tximeta object.")
  }

  # Shuffle only C1/H1 membership. H3 and H24 labels remain unchanged so each
  # fit uses the same complete time-course model as the manuscript analysis.
  shuffled_condition <- as.character(input_coldata$condition)
  names(shuffled_condition) <- rownames(input_coldata)
  shuffled_condition[split$A] <- "C1"
  shuffled_condition[split$B] <- "H1"

  original_levels <- levels(factor(input_coldata$condition))
  preferred_levels <- c("C1", "H1", "H3", "H24")
  condition_levels <- c(
    intersect(preferred_levels, unique(shuffled_condition)),
    setdiff(original_levels, preferred_levels)
  )
  input_coldata$condition <- stats::relevel(
    factor(shuffled_condition, levels = unique(condition_levels)),
    ref = "C1"
  )

  # Start from the same pre-DESeq, gene-summarized tximeta object used by
  # run_deseq.R. This retains tximeta's counts and transcript-length offsets.
  input_se <- gene_se
  SummarizedExperiment::colData(input_se) <- S4Vectors::DataFrame(
    input_coldata
  )
  dds <- DESeq2::DESeqDataSet(
    input_se,
    design = DESeq2::design(source_dds)
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- DESeq2::results(
    dds,
    contrast = c("condition", "H1", "C1"),
    alpha = alpha
  )

  tested <- !is.na(res$padj)
  sig <- tested & res$padj < alpha
  tested_ids <- unique(.sup_strip_version(rownames(res)[tested]))
  sig_ids <- .sup_strip_version(rownames(res)[sig])

  list(
    summary = .sup_summary_row(
      split_id,
      "DESeq2",
      "gene",
      length(tested_ids),
      length(unique(sig_ids)),
      sprintf("BH FDR < %.3g", alpha)
    ),
    sig_genes = list(DESeq2 = unique(sig_ids)),
    meta = list(
      analysis_scope = "Complete time course; only C1/H1 labels shuffled",
      design = deparse(DESeq2::design(source_dds)),
      input = "Pre-DESeq gene-level tximeta SummarizedExperiment",
      retained_assays = SummarizedExperiment::assayNames(input_se)
    )
  )
}

.sup_run_dexseq_split <- function(
  split_id,
  split,
  dexseq_results,
  alpha
) {
  .sup_check_pkg(c("DEXSeq", "SummarizedExperiment", "S4Vectors"))
  source_se <- dexseq_results$by_tp$H1$se
  if (is.null(source_se)) {
    stop(
      "Cached DEXSeq H1 result has no $se object. Use a summarizeOverlaps ",
      "DEXSeq cache or rerun DEXSeq once with that backend."
    )
  }

  sample_ids <- c(split$A, split$B)
  if (!all(sample_ids %in% colnames(source_se))) {
    stop("DEXSeq cached exon-bin SE is missing permutation samples.")
  }
  se <- source_se[, sample_ids, drop = FALSE]
  dex_sample_data <- data.frame(
    row.names = sample_ids,
    sample = factor(sample_ids),
    condition = factor(
      rep(c("A", "B"), c(length(split$A), length(split$B))),
      levels = c("A", "B")
    ),
    libType = factor("paired-end")
  )
  SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(dex_sample_data)

  dxd <- suppressMessages(suppressWarnings(
    DEXSeq::DEXSeqDataSetFromSE(
      se,
      design = ~sample + exon + condition:exon
    )
  ))
  dxd <- DEXSeq::estimateSizeFactors(dxd)
  dxd <- DEXSeq::estimateDispersions(dxd)
  dxd <- DEXSeq::testForDEU(dxd)
  dxr <- DEXSeq::DEXSeqResults(dxd)

  exon_padj <- dxr$padj
  exon_tested <- !is.na(exon_padj)
  exon_sig <- exon_tested & exon_padj <= alpha

  gene_q <- DEXSeq::perGeneQValue(dxr)
  gene_tested <- !is.na(gene_q)
  gene_sig <- gene_tested & gene_q <= alpha

  list(
    summary = rbind(
      .sup_summary_row(
        split_id,
        "DEXSeq",
        "gene",
        sum(gene_tested),
        sum(gene_sig),
        sprintf("perGeneQValue <= %.3g", alpha)
      ),
      .sup_summary_row(
        split_id,
        "DEXSeq",
        "exon_bin",
        sum(exon_tested),
        sum(exon_sig),
        sprintf("BH FDR <= %.3g", alpha)
      )
    ),
    sig_genes = list(
      DEXSeq = unique(.sup_strip_version(names(gene_q)[gene_sig]))
    )
  )
}

.sup_run_drimseq_split <- function(
  split_id,
  split,
  drimseq_results,
  alpha,
  BPPARAM,
  seed,
  precision_mode = c("reuse_observed", "refit")
) {
  precision_mode <- match.arg(precision_mode)
  .sup_check_pkg(c("DRIMSeq", "stageR"))
  counts_df <- drimseq_results$by_tp$H1$counts
  if (is.null(counts_df)) {
    stop("Expected drimseq_results$by_tp$H1$counts.")
  }

  sample_ids <- c(split$A, split$B)
  if (!all(sample_ids %in% colnames(counts_df))) {
    stop("DRIMSeq cached count table is missing permutation samples.")
  }
  counts_df <- counts_df[, c("gene_id", "feature_id", sample_ids), drop = FALSE]
  keep <- rowSums(counts_df[, sample_ids, drop = FALSE]) > 0
  counts_df <- counts_df[keep, , drop = FALSE]

  samples_df <- data.frame(
    sample_id = sample_ids,
    condition = factor(
      rep(c("A", "B"), c(length(split$A), length(split$B))),
      levels = c("A", "B")
    ),
    stringsAsFactors = FALSE
  )
  d <- DRIMSeq::dmDSdata(counts = counts_df, samples = samples_df)
  d <- DRIMSeq::dmFilter(
    d,
    min_samps_feature_expr = 3,
    min_samps_gene_expr = 6,
    min_feature_expr = 10,
    min_gene_expr = 10
  )
  design <- stats::model.matrix(~condition, data = DRIMSeq::samples(d))

  if (precision_mode == "reuse_observed") {
    observed_dm <- drimseq_results$by_tp$H1$dm
    if (is.null(observed_dm)) {
      stop(
        "Cannot reuse DRIMSeq precision: ",
        "drimseq_results$by_tp$H1$dm is missing."
      )
    }

    current_counts <- DRIMSeq::counts(d)
    observed_counts <- DRIMSeq::counts(observed_dm)
    current_keys <- paste(current_counts$gene_id, current_counts$feature_id, sep = "\r")
    observed_keys <- paste(
      observed_counts$gene_id,
      observed_counts$feature_id,
      sep = "\r"
    )
    if (!setequal(current_keys, observed_keys)) {
      stop(
        "Cannot reuse DRIMSeq precision because the filtered gene/transcript ",
        "universe differs from the cached observed H1 analysis."
      )
    }

    observed_genewise <- DRIMSeq::genewise_precision(observed_dm)
    current_gene_ids <- unique(current_counts$gene_id)
    missing_precision <- setdiff(
      current_gene_ids,
      observed_genewise$gene_id
    )
    if (length(missing_precision) > 0) {
      stop(
        "Cached observed DRIMSeq object lacks gene-wise precision for ",
        length(missing_precision), " filtered genes."
      )
    }

    # Construct a precision-class object for the relabeled samples without
    # re-running the expensive common and gene-wise precision optimizations.
    # Mean expression is count-derived and inexpensive to recalculate.
    d <- DRIMSeq::dmPrecision(
      d,
      design = design,
      mean_expression = TRUE,
      common_precision = FALSE,
      genewise_precision = FALSE,
      BPPARAM = BiocParallel::SerialParam()
    )
    observed_common <- DRIMSeq::common_precision(observed_dm)
    if (length(observed_common) == 1L) {
      DRIMSeq::common_precision(d) <- observed_common
    }
    DRIMSeq::genewise_precision(d) <- observed_genewise
    message(
      "[DRIMSeq] ", split_id,
      ": reused observed H1 gene-wise precision estimates for ",
      length(current_gene_ids), " genes"
    )
  } else {
    message(
      "[DRIMSeq] ", split_id,
      ": estimating common and gene-wise precision"
    )
    set.seed(seed)
    d <- DRIMSeq::dmPrecision(d, design = design, BPPARAM = BPPARAM)
  }

  message("[DRIMSeq] ", split_id, ": fitting relabeled models")
  d <- DRIMSeq::dmFit(d, design = design, BPPARAM = BPPARAM)
  message("[DRIMSeq] ", split_id, ": testing differential transcript usage")
  d <- DRIMSeq::dmTest(d, coef = "conditionB")

  gene_res <- as.data.frame(DRIMSeq::results(d))
  tx_res <- as.data.frame(DRIMSeq::results(d, level = "feature"))

  p_screen <- gene_res$pvalue
  p_screen[is.na(p_screen)] <- 1
  names(p_screen) <- gene_res$gene_id
  p_confirmation <- matrix(tx_res$pvalue, ncol = 1)
  p_confirmation[is.na(p_confirmation)] <- 1
  rownames(p_confirmation) <- tx_res$feature_id

  stage <- stageR::stageRTx(
    pScreen = p_screen,
    pConfirmation = p_confirmation,
    pScreenAdjusted = FALSE,
    tx2gene = tx_res[, c("feature_id", "gene_id")]
  )
  stage <- stageR::stageWiseAdjustment(stage, method = "dtu", alpha = alpha)
  sig_gene_ids <- suppressMessages(
    rownames(as.data.frame(stageR::getSignificantGenes(stage)))
  )
  sig_tx_ids <- suppressMessages(
    rownames(as.data.frame(stageR::getSignificantTx(stage)))
  )

  list(
    summary = rbind(
      .sup_summary_row(
        split_id,
        "DRIMSeq",
        "gene",
        sum(!is.na(gene_res$pvalue)),
        length(unique(sig_gene_ids)),
        sprintf("stageR OFDR %.3g", alpha)
      ),
      .sup_summary_row(
        split_id,
        "DRIMSeq",
        "transcript",
        sum(!is.na(tx_res$pvalue)),
        length(unique(sig_tx_ids)),
        sprintf("stageR OFDR %.3g", alpha)
      )
    ),
    sig_genes = list(
      DRIMSeq = unique(.sup_strip_version(sig_gene_ids))
    )
  )
}

.sup_read_suppa_isotpm <- function(path) {
  if (!file.exists(path)) stop("SUPPA2 master isoTPM file not found: ", path)
  con <- file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  sample_names <- strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]]
  x <- utils::read.delim(
    con,
    header = FALSE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (ncol(x) != length(sample_names) + 1) {
    stop("Unexpected SUPPA2 master isoTPM dimensions: ", path)
  }
  tx_ids <- x[[1]]
  mat <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- tx_ids
  colnames(mat) <- sample_names
  mat
}

.sup_write_suppa_isotpm <- function(mat, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(paste(colnames(mat), collapse = "\t"), con = con)
  utils::write.table(
    data.frame(transcript_id = rownames(mat), mat, check.names = FALSE),
    file = con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  invisible(path)
}

.sup_find_suppa_ioe <- function(ioe_file = NULL) {
  if (!is.null(ioe_file)) return(ioe_file)
  candidates <- list.files(
    file.path(PROJECT_ROOT, "resources", "suppa_events"),
    pattern = "\\.events\\.ioe$",
    full.names = TRUE
  )
  if (length(candidates) != 1) {
    stop(
      "Could not uniquely resolve the SUPPA2 master IOE file. ",
      "Pass suppa_ioe_file explicitly."
    )
  }
  candidates[1]
}

.sup_system2 <- function(command, args, log_file) {
  status <- system2(
    command = command,
    args = vapply(args, shQuote, character(1)),
    stdout = log_file,
    stderr = log_file
  )
  if (length(status) != 1 || is.na(status) || status != 0) {
    stop(
      "Command failed with status ", status, ": ", command,
      "\nSee log: ", log_file
    )
  }
  invisible(TRUE)
}

.sup_resolve_suppa_command <- function(
  suppa_executable,
  suppa_runner = NULL
) {
  if (!is.null(suppa_runner) && nzchar(suppa_runner)) {
    runner <- Sys.which(suppa_runner)
    if (!nzchar(runner) && file.exists(suppa_runner)) {
      runner <- normalizePath(suppa_runner)
    }
    if (!nzchar(runner)) {
      stop("SUPPA2 environment runner not found: ", suppa_runner)
    }
    if (file.access(runner, mode = 1) != 0) {
      bash <- Sys.which("bash")
      if (!nzchar(bash)) {
        stop(
          "SUPPA2 environment runner is not executable and bash was not ",
          "found: ", runner
        )
      }
      return(list(
        command = bash,
        args_prefix = c(runner, suppa_executable),
        description = paste(bash, runner, suppa_executable)
      ))
    }
    return(list(
      command = runner,
      args_prefix = suppa_executable,
      description = paste(runner, suppa_executable)
    ))
  }

  executable <- Sys.which(suppa_executable)
  if (!nzchar(executable) && file.exists(suppa_executable)) {
    executable <- normalizePath(suppa_executable)
  }
  if (!nzchar(executable)) {
    stop(
      "SUPPA2 executable not found on the current PATH: ",
      suppa_executable,
      "\nEither pass suppa_runner=file.path(PROJECT_ROOT, 'scripts', 'run') ",
      "to use the dedicated conda environment, or supply the full path to ",
      "suppa.py."
    )
  }
  list(
    command = executable,
    args_prefix = character(0),
    description = executable
  )
}

.sup_read_dpsi <- function(path) {
  x <- utils::read.delim(
    path,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (ncol(x) < 2) stop("Unexpected SUPPA2 dpsi file: ", path)
  cn <- colnames(x)
  p_col <- grep("p[-._]?val|pvalue", cn, ignore.case = TRUE, value = TRUE)
  if (length(p_col) == 0) p_col <- cn[2]

  data.frame(
    event_id = rownames(x),
    gene_id = .sup_strip_version(sub(";.*$", "", rownames(x))),
    pvalue = suppressWarnings(as.numeric(x[[p_col[1]]])),
    stringsAsFactors = FALSE
  )
}

.sup_run_suppa_split <- function(
  split_id,
  split,
  alpha,
  master_tpm,
  ioe_file,
  suppa_executable,
  suppa_runner,
  split_dir
) {
  suppa_command <- .sup_resolve_suppa_command(
    suppa_executable = suppa_executable,
    suppa_runner = suppa_runner
  )
  if (!file.exists(ioe_file)) stop("SUPPA2 IOE file not found: ", ioe_file)

  sample_ids <- c(split$A, split$B)
  if (!all(sample_ids %in% colnames(master_tpm))) {
    stop("SUPPA2 master isoTPM matrix is missing permutation samples.")
  }
  dir.create(split_dir, recursive = TRUE, showWarnings = FALSE)

  a_tpm <- file.path(split_dir, "group_A_iso_tpm.txt")
  b_tpm <- file.path(split_dir, "group_B_iso_tpm.txt")
  a_prefix <- file.path(split_dir, "group_A")
  b_prefix <- file.path(split_dir, "group_B")
  a_psi <- paste0(a_prefix, ".psi")
  b_psi <- paste0(b_prefix, ".psi")
  out_prefix <- file.path(split_dir, "diffSplice")
  dpsi_file <- paste0(out_prefix, ".dpsi")

  .sup_write_suppa_isotpm(master_tpm[, split$A, drop = FALSE], a_tpm)
  .sup_write_suppa_isotpm(master_tpm[, split$B, drop = FALSE], b_tpm)

  .sup_system2(
    suppa_command$command,
    c(
      suppa_command$args_prefix,
      "psiPerEvent", "-i", ioe_file, "-e", a_tpm, "-o", a_prefix
    ),
    file.path(split_dir, "suppa_group_A.log")
  )
  .sup_system2(
    suppa_command$command,
    c(
      suppa_command$args_prefix,
      "psiPerEvent", "-i", ioe_file, "-e", b_tpm, "-o", b_prefix
    ),
    file.path(split_dir, "suppa_group_B.log")
  )
  if (!file.exists(a_psi) || !file.exists(b_psi)) {
    stop("SUPPA2 did not create both PSI files for ", split_id)
  }

  .sup_system2(
    suppa_command$command,
    c(
      suppa_command$args_prefix,
      "diffSplice",
      "-m", "empirical",
      "-gc",
      "-s",
      "-i", ioe_file,
      "-p", a_psi, b_psi,
      "-e", a_tpm, b_tpm,
      "-o", out_prefix
    ),
    file.path(split_dir, "suppa_diffSplice.log")
  )
  if (!file.exists(dpsi_file)) {
    stop("SUPPA2 dpsi output not found for ", split_id)
  }

  res <- .sup_read_dpsi(dpsi_file)
  res$padj_BH <- stats::p.adjust(res$pvalue, method = "BH")
  tested <- !is.na(res$pvalue)
  raw_sig <- tested & res$pvalue < alpha
  bh_sig <- tested & res$padj_BH < alpha

  list(
    summary = rbind(
      .sup_summary_row(
        split_id,
        "SUPPA2_rawP",
        "gene",
        length(unique(res$gene_id[tested])),
        length(unique(res$gene_id[raw_sig])),
        sprintf("any event unadjusted p < %.3g", alpha)
      ),
      .sup_summary_row(
        split_id,
        "SUPPA2_rawP",
        "event",
        sum(tested),
        sum(raw_sig),
        sprintf("unadjusted p < %.3g", alpha)
      ),
      .sup_summary_row(
        split_id,
        "SUPPA2_BH",
        "gene",
        length(unique(res$gene_id[tested])),
        length(unique(res$gene_id[bh_sig])),
        sprintf("any event BH FDR < %.3g", alpha)
      ),
      .sup_summary_row(
        split_id,
        "SUPPA2_BH",
        "event",
        sum(tested),
        sum(bh_sig),
        sprintf("BH FDR < %.3g", alpha)
      )
    ),
    sig_genes = list(
      SUPPA2_rawP = unique(res$gene_id[raw_sig]),
      SUPPA2_BH = unique(res$gene_id[bh_sig])
    )
  )
}

.sup_observed_results <- function(
  deseq_results,
  dexseq_results,
  drimseq_results,
  suppa_results,
  alpha
) {
  rows <- list()
  sig <- list()

  if (!is.null(deseq_results)) {
    full <- deseq_results$results$gene_full_by_tp$H1
    ids <- unique(deseq_results$results$gene_sig_by_tp$H1$ensgene)
    rows[[length(rows) + 1L]] <- .sup_summary_row(
      "observed",
      "DESeq2",
      "gene",
      sum(!is.na(full$padj)),
      length(ids),
      sprintf("BH FDR < %.3g", alpha)
    )
    sig$DESeq2 <- .sup_strip_version(ids)
  }

  if (!is.null(dexseq_results)) {
    gf <- dexseq_results$results$gene_full_by_tp$H1
    ef <- dexseq_results$results$exon_full_by_tp$H1
    gids <- unique(dexseq_results$results$gene_sig_by_tp$H1$ensgene)
    rows[[length(rows) + 1L]] <- rbind(
      .sup_summary_row(
        "observed",
        "DEXSeq",
        "gene",
        sum(!is.na(gf$padj)),
        length(gids),
        sprintf("perGeneQValue <= %.3g", alpha)
      ),
      .sup_summary_row(
        "observed",
        "DEXSeq",
        "exon_bin",
        sum(!is.na(ef$padj)),
        nrow(dexseq_results$results$exon_sig_by_tp$H1),
        sprintf("BH FDR <= %.3g", alpha)
      )
    )
    sig$DEXSeq <- .sup_strip_version(gids)
  }

  if (!is.null(drimseq_results)) {
    gf <- drimseq_results$results$gene_full_by_tp$H1
    tf <- drimseq_results$results$tx_full_by_tp$H1
    gids <- unique(drimseq_results$results$gene_sig_by_tp$H1$ensgene)
    tids <- unique(drimseq_results$results$tx_sig_by_tp$H1$feature_id)
    rows[[length(rows) + 1L]] <- rbind(
      .sup_summary_row(
        "observed",
        "DRIMSeq",
        "gene",
        sum(!is.na(gf$pvalue)),
        length(gids),
        sprintf("stageR OFDR %.3g", alpha)
      ),
      .sup_summary_row(
        "observed",
        "DRIMSeq",
        "transcript",
        sum(!is.na(tf$pvalue)),
        length(tids),
        sprintf("stageR OFDR %.3g", alpha)
      )
    )
    sig$DRIMSeq <- .sup_strip_version(gids)
  }

  if (!is.null(suppa_results)) {
    full <- suppa_results$results$event_full_by_tp$H1
    pvalue <- full$pvalue
    padj <- stats::p.adjust(pvalue, method = "BH")
    tested <- !is.na(pvalue)
    raw_sig <- tested & pvalue < alpha
    bh_sig <- tested & padj < alpha
    gene_ids <- .sup_strip_version(full$ensgene)

    rows[[length(rows) + 1L]] <- rbind(
      .sup_summary_row(
        "observed",
        "SUPPA2_rawP",
        "gene",
        length(unique(gene_ids[tested])),
        length(unique(gene_ids[raw_sig])),
        sprintf("any event unadjusted p < %.3g", alpha)
      ),
      .sup_summary_row(
        "observed",
        "SUPPA2_rawP",
        "event",
        sum(tested),
        sum(raw_sig),
        sprintf("unadjusted p < %.3g", alpha)
      ),
      .sup_summary_row(
        "observed",
        "SUPPA2_BH",
        "gene",
        length(unique(gene_ids[tested])),
        length(unique(gene_ids[bh_sig])),
        sprintf("any event BH FDR < %.3g", alpha)
      ),
      .sup_summary_row(
        "observed",
        "SUPPA2_BH",
        "event",
        sum(tested),
        sum(bh_sig),
        sprintf("BH FDR < %.3g", alpha)
      )
    )
    sig$SUPPA2_rawP <- unique(gene_ids[raw_sig])
    sig$SUPPA2_BH <- unique(gene_ids[bh_sig])
  }

  list(summary = do.call(rbind, rows), sig_genes = sig)
}

.sup_error_result <- function(split_id, tool, message_text) {
  levels <- switch(
    tool,
    DESeq2 = "gene",
    DEXSeq = c("gene", "exon_bin"),
    DRIMSeq = c("gene", "transcript"),
    SUPPA2 = c("gene", "event", "gene", "event")
  )
  tool_names <- if (tool == "SUPPA2") {
    c("SUPPA2_rawP", "SUPPA2_rawP", "SUPPA2_BH", "SUPPA2_BH")
  } else {
    rep(tool, length(levels))
  }

  rows <- do.call(rbind, lapply(seq_along(levels), function(i) {
    .sup_summary_row(
      split_id,
      tool_names[i],
      levels[i],
      NA_integer_,
      NA_integer_,
      "error",
      status = message_text
    )
  }))
  list(summary = rows, sig_genes = list())
}

.sup_plot_metrics <- function(summary_df, suppa_threshold) {
  suppa_key <- if (suppa_threshold == "paper_raw_p") {
    "SUPPA2_rawP"
  } else {
    "SUPPA2_BH"
  }

  get_row <- function(tool, label) {
    z <- summary_df[
      summary_df$tool == tool &
        summary_df$level == "gene" &
        summary_df$status == "ok",
      ,
      drop = FALSE
    ]
    data.frame(
      metric = label,
      n_tested = if (nrow(z) == 1) z$n_tested else NA_integer_,
      n_significant = if (nrow(z) == 1) z$n_significant else NA_integer_,
      stringsAsFactors = FALSE
    )
  }

  rbind(
    get_row("DESeq2", "DESeq2"),
    get_row("DEXSeq", "DEXSeq"),
    get_row("DRIMSeq", "DRIMSeq"),
    get_row(suppa_key, "SUPPA2")
  )
}

.sup_run_permutations <- function(
  deseq_results,
  gene_se,
  dexseq_results,
  drimseq_results,
  suppa_results,
  permutation_tools = c("DESeq2", "DEXSeq", "DRIMSeq", "SUPPA2"),
  n_null_splits = Inf,
  alpha = 0.10,
  suppa_threshold = c("paper_raw_p", "BH"),
  suppa_master_tpm_file = NULL,
  suppa_ioe_file = NULL,
  suppa_executable = "suppa.py",
  suppa_runner = file.path(PROJECT_ROOT, "scripts", "run"),
  drimseq_precision_mode = c("reuse_observed", "refit"),
  workers = 1L,
  seed = 333L,
  cache_dir,
  force = FALSE
) {
  suppa_threshold <- match.arg(suppa_threshold)
  drimseq_precision_mode <- match.arg(drimseq_precision_mode)
  valid_tools <- c("DESeq2", "DEXSeq", "DRIMSeq", "SUPPA2")
  unknown_tools <- setdiff(permutation_tools, valid_tools)
  if (length(unknown_tools) > 0) {
    stop("Unknown permutation tool(s): ", paste(unknown_tools, collapse = ", "))
  }

  split_info <- .sup_make_balanced_splits(deseq_results)
  null_splits <- .sup_select_null_splits(
    split_info$null,
    n_null_splits = n_null_splits,
    seed = seed
  )
  membership <- .sup_membership_table(
    c(list(observed = split_info$observed), null_splits),
    split_info$coldata
  )
  split_descriptors <- .sup_split_descriptors(membership)
  descriptor_idx <- match(membership$split_id, split_descriptors$split_id)
  membership$split_short <- split_descriptors$split_short[descriptor_idx]
  membership$swap_label <- split_descriptors$swap_label[descriptor_idx]

  .sup_check_pkg("BiocParallel")
  workers <- max(1L, as.integer(workers))
  BPPARAM <- if (.Platform$OS.type == "windows") {
    BiocParallel::SnowParam(workers = workers, progressbar = FALSE)
  } else {
    BiocParallel::MulticoreParam(workers = workers, progressbar = FALSE)
  }

  master_tpm <- NULL
  ioe_file <- NULL
  if ("SUPPA2" %in% permutation_tools) {
    suppa_master_tpm_file <- suppa_master_tpm_file %||% file.path(
      .get_results_dir(),
      "analysis", "suppa", "standard", "C1", "iso_tpm", "iso_tpm.txt"
    )
    master_tpm <- .sup_read_suppa_isotpm(suppa_master_tpm_file)
    ioe_file <- .sup_find_suppa_ioe(suppa_ioe_file)
  }

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  observed <- .sup_observed_results(
    deseq_results = if ("DESeq2" %in% permutation_tools) deseq_results else NULL,
    dexseq_results = if ("DEXSeq" %in% permutation_tools) dexseq_results else NULL,
    drimseq_results = if ("DRIMSeq" %in% permutation_tools) drimseq_results else NULL,
    suppa_results = if ("SUPPA2" %in% permutation_tools) suppa_results else NULL,
    alpha = alpha
  )

  deseq_observed_verification <- NULL
  if ("DESeq2" %in% permutation_tools) {
    observed_deseq <- .sup_run_deseq_split(
      split_id = "observed",
      split = split_info$observed,
      deseq_results = deseq_results,
      gene_se = gene_se,
      alpha = alpha
    )
    observed$summary <- rbind(
      observed$summary[observed$summary$tool != "DESeq2", , drop = FALSE],
      observed_deseq$summary
    )
    observed$sig_genes$DESeq2 <- observed_deseq$sig_genes$DESeq2

    original_n <- nrow(
      deseq_results$results$gene_sig_by_tp$H1
    )
    refit_n <- observed_deseq$summary$n_significant[
      observed_deseq$summary$tool == "DESeq2" &
        observed_deseq$summary$level == "gene"
    ]
    deseq_observed_verification <- data.frame(
      original_cached_deg_count = as.integer(original_n),
      observed_refit_deg_count = as.integer(refit_n),
      exact_match = length(refit_n) == 1 && original_n == refit_n,
      stringsAsFactors = FALSE
    )
    .sup_write_tsv(
      deseq_observed_verification,
      file.path(cache_dir, "deseq_observed_verification.tsv")
    )
    if (!isTRUE(deseq_observed_verification$exact_match)) {
      stop(
        "Observed DESeq2 permutation refit returned ", refit_n,
        " DEGs, but the cached manuscript result contains ", original_n,
        ". The permutation run has been stopped before fitting null splits. ",
        "See: ", file.path(cache_dir, "deseq_observed_verification.tsv")
      )
    } else {
      message(
        "[DESeq2 verification] observed refit reproduced the cached result: ",
        original_n, " significant genes."
      )
    }
  }

  all_summary <- list(observed$summary)
  all_plot_metrics <- list(
    transform(
      .sup_plot_metrics(
        observed$summary,
        suppa_threshold
      ),
      split_id = "observed",
      split_type = "Observed"
    )
  )
  all_sig <- list(observed = observed$sig_genes)

  for (split_index in seq_along(null_splits)) {
    split_id <- names(null_splits)[split_index]
    split <- null_splits[[split_index]]
    message(
      "\n[Permutation] ", split_id,
      " (", split_index, "/", length(null_splits), ")"
    )

    split_results <- list()
    for (tool in permutation_tools) {
      tool_cache <- file.path(
        cache_dir,
        split_id,
        paste0(
          tolower(tool), "_",
          unname(.SUP_PERMUTATION_CACHE_VERSION[tool]),
          if (tool == "DRIMSeq") {
            paste0("_", drimseq_precision_mode)
          } else {
            ""
          },
          ".rds"
        )
      )
      dir.create(dirname(tool_cache), recursive = TRUE, showWarnings = FALSE)

      if (file.exists(tool_cache) && !force) {
        cached_tool <- readRDS(tool_cache)
        cached_ok <- !is.null(cached_tool$summary) &&
          all(cached_tool$summary$status == "ok")
        if (cached_ok) {
          message("[", tool, "] loading cached split: ", split_id)
          split_results[[tool]] <- cached_tool
          next
        }
        message(
          "[", tool, "] retrying a previously failed cached split: ", split_id
        )
      }

      message("[", tool, "] fitting split: ", split_id)
      result <- tryCatch(
        switch(
          tool,
          DESeq2 = .sup_run_deseq_split(
            split_id, split, deseq_results, gene_se, alpha
          ),
          DEXSeq = .sup_run_dexseq_split(
            split_id, split, dexseq_results, alpha
          ),
          DRIMSeq = .sup_run_drimseq_split(
            split_id,
            split,
            drimseq_results,
            alpha,
            BPPARAM,
            seed + split_index,
            precision_mode = drimseq_precision_mode
          ),
          SUPPA2 = .sup_run_suppa_split(
            split_id = split_id,
            split = split,
            alpha = alpha,
            master_tpm = master_tpm,
            ioe_file = ioe_file,
            suppa_executable = suppa_executable,
            suppa_runner = suppa_runner,
            split_dir = file.path(cache_dir, split_id, "suppa2_files")
          )
        ),
        error = function(e) {
          warning("[", tool, "] ", split_id, " failed: ", conditionMessage(e))
          .sup_error_result(split_id, tool, conditionMessage(e))
        }
      )
      saveRDS(result, tool_cache)
      split_results[[tool]] <- result
    }

    split_summary <- do.call(
      rbind,
      lapply(split_results, `[[`, "summary")
    )
    split_sig <- unlist(
      lapply(split_results, `[[`, "sig_genes"),
      recursive = FALSE
    )

    all_summary[[length(all_summary) + 1L]] <- split_summary
    all_plot_metrics[[length(all_plot_metrics) + 1L]] <- transform(
      .sup_plot_metrics(split_summary, suppa_threshold),
      split_id = split_id,
      split_type = "Permuted"
    )
    all_sig[[split_id]] <- split_sig
  }

  summary_df <- do.call(rbind, all_summary)
  plot_metrics <- do.call(rbind, all_plot_metrics)
  rownames(summary_df) <- NULL
  rownames(plot_metrics) <- NULL
  descriptor_idx <- match(plot_metrics$split_id, split_descriptors$split_id)
  plot_metrics$split_short <- split_descriptors$split_short[descriptor_idx]
  plot_metrics$swap_label <- split_descriptors$swap_label[descriptor_idx]
  plot_metrics$swapped_C1_sample <-
    split_descriptors$swapped_C1_sample[descriptor_idx]
  plot_metrics$swapped_H1_sample <-
    split_descriptors$swapped_H1_sample[descriptor_idx]

  expected_null <- length(split_info$null)
  randomization <- do.call(rbind, lapply(unique(plot_metrics$metric), function(metric) {
    z <- plot_metrics[plot_metrics$metric == metric, , drop = FALSE]
    obs <- z$n_significant[z$split_type == "Observed"]
    null <- z$n_significant[z$split_type == "Permuted"]
    complete_exhaustive <- length(null) == expected_null &&
      length(obs) == 1 &&
      !is.na(obs) &&
      all(!is.na(null))

    data.frame(
      metric = metric,
      observed = if (length(obs) == 1) obs else NA_integer_,
      n_null_splits = sum(!is.na(null)),
      null_min = if (any(!is.na(null))) min(null, na.rm = TRUE) else NA_real_,
      null_median = if (any(!is.na(null))) stats::median(null, na.rm = TRUE) else NA_real_,
      null_max = if (any(!is.na(null))) max(null, na.rm = TRUE) else NA_real_,
      n_null_ge_observed = if (length(obs) == 1 && any(!is.na(null))) {
        sum(null >= obs, na.rm = TRUE)
      } else {
        NA_integer_
      },
      observed_to_null_median_ratio = if (
        length(obs) == 1 && any(!is.na(null))
      ) {
        null_median <- stats::median(null, na.rm = TRUE)
        if (null_median == 0 && obs > 0) Inf else obs / null_median
      } else {
        NA_real_
      },
      upper_tail_partition_fraction = if (complete_exhaustive) {
        (1 + sum(null >= obs)) / (1 + length(null))
      } else {
        NA_real_
      },
      exact_upper_tail_p = if (complete_exhaustive) {
        (1 + sum(null >= obs)) / (1 + length(null))
      } else {
        NA_real_
      },
      exhaustive = complete_exhaustive,
      stringsAsFactors = FALSE
    )
  }))

  if (
    !is.null(deseq_observed_verification) &&
      !isTRUE(deseq_observed_verification$exact_match)
  ) {
    idx <- randomization$metric == "DESeq2"
    randomization$exact_upper_tail_p[idx] <- NA_real_
    randomization$exhaustive[idx] <- FALSE
  }

  if (
    "DRIMSeq" %in% permutation_tools &&
      drimseq_precision_mode == "reuse_observed"
  ) {
    idx <- randomization$metric == "DRIMSeq"
    randomization$exact_upper_tail_p[idx] <- NA_real_
    randomization$exhaustive[idx] <- FALSE
  }

  list(
    meta = list(
      analysis_version = .SUP_FIG_X_ANALYSIS_VERSION,
      alpha = alpha,
      seed = seed,
      workers = workers,
      permutation_tools = permutation_tools,
      suppa_threshold = suppa_threshold,
      suppa_runner = suppa_runner,
      suppa_executable = suppa_executable,
      drimseq_precision_mode = drimseq_precision_mode,
      drimseq_precision_note = if (
        drimseq_precision_mode == "reuse_observed"
      ) {
        paste(
          "Gene-wise precision estimates were fixed to the cached observed",
          "C1-versus-H1 fit. Null counts are a conditional/descriptive",
          "sensitivity analysis; no exact DRIMSeq permutation p-value is reported."
        )
      } else {
        paste(
          "Common and gene-wise precision were re-estimated for every split;",
          "the DRIMSeq randomization p-value is exhaustive when all splits succeed."
        )
      },
      n_unique_balanced_partitions = 10L,
      n_alternative_null_partitions = expected_null,
      n_null_splits_run = length(null_splits),
      exhaustive = length(null_splits) == expected_null,
      deseq_model = paste(
        "Complete original time course; H3/H24 fixed;",
        "C1/H1 membership permuted among the six early samples."
      ),
      interpretation = paste(
        "Internal label-permutation diagnostic using the six C1/H1 samples;",
        "not an independent biological negative-control experiment."
      )
    ),
    membership = membership,
    split_descriptors = split_descriptors,
    tool_summary = summary_df,
    plot_metrics = plot_metrics,
    randomization_summary = randomization,
    deseq_observed_verification = deseq_observed_verification,
    significant_genes = all_sig
  )
}

# --------------------------------------------------
# Main analysis wrapper
# --------------------------------------------------

run_sup_fig_x_analysis <- function(
  deseq_results = NULL,
  gene_se = NULL,
  dexseq_results = NULL,
  drimseq_results = NULL,
  suppa_results = NULL,
  transcript_se = NULL,
  tx2gene = NULL,
  outdir = NULL,
  pca_feature_numbers = c(500L, 1000L, 5000L),
  gene_pca_display = "all_tested",
  transcript_pca_display = "all_expressed",
  transcript_min_tpm = 1,
  transcript_min_samples = 3L,
  transcript_pseudocount = 0.5,
  include_transcript_gsea = TRUE,
  hypoxia_genes = NULL,
  run_pairwise_psi = TRUE,
  pairwise_psi_cutoff = 0.10,
  suppa_psi_dir = NULL,
  run_permutations = TRUE,
  permutation_tools = c("DESeq2", "DEXSeq", "DRIMSeq", "SUPPA2"),
  n_null_splits = Inf,
  padj_cutoff = 0.10,
  suppa_threshold = c("paper_raw_p", "BH"),
  suppa_master_tpm_file = NULL,
  suppa_ioe_file = NULL,
  suppa_executable = "suppa.py",
  suppa_runner = file.path(PROJECT_ROOT, "scripts", "run"),
  drimseq_precision_mode = c("reuse_observed", "refit"),
  workers = 1L,
  seed = 333L,
  save_tables = TRUE,
  force = FALSE
) {
  suppa_threshold <- match.arg(suppa_threshold)
  drimseq_precision_mode <- match.arg(drimseq_precision_mode)
  if (is.null(outdir)) {
    outdir <- file.path(.get_results_dir(), "plots", "sup_fig_x")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  cache_file <- file.path(outdir, "sup_fig_x_analysis.rds")

  if (file.exists(cache_file) && !force) {
    cached <- readRDS(cache_file)
    cache_matches <- identical(
      cached$meta$analysis_version,
      .SUP_FIG_X_ANALYSIS_VERSION
    ) &&
      identical(cached$meta$pca_feature_numbers, pca_feature_numbers) &&
      identical(cached$meta$gene_pca_display, gene_pca_display) &&
      identical(
        cached$meta$transcript_pca_display,
        transcript_pca_display
      ) &&
      identical(cached$meta$transcript_min_tpm, transcript_min_tpm) &&
      identical(
        cached$meta$transcript_min_samples,
        transcript_min_samples
      ) &&
      identical(
        cached$meta$transcript_pseudocount,
        transcript_pseudocount
      ) &&
      identical(
        cached$meta$include_transcript_gsea,
        include_transcript_gsea
      ) &&
      identical(cached$meta$run_pairwise_psi, run_pairwise_psi) &&
      identical(cached$meta$pairwise_psi_cutoff, pairwise_psi_cutoff) &&
      identical(cached$meta$suppa_psi_dir, suppa_psi_dir) &&
      identical(cached$meta$run_permutations, run_permutations) &&
      identical(cached$meta$padj_cutoff, padj_cutoff) &&
      identical(cached$meta$permutation_tools, permutation_tools) &&
      identical(cached$meta$n_null_splits, n_null_splits) &&
      identical(cached$meta$suppa_threshold, suppa_threshold) &&
      identical(cached$meta$suppa_executable, suppa_executable) &&
      identical(cached$meta$suppa_runner, suppa_runner) &&
      identical(
        cached$meta$drimseq_precision_mode,
        drimseq_precision_mode
      )

    if (isTRUE(cache_matches)) {
      message("Loading cached supplemental analysis: ", cache_file)
      return(cached)
    }
    message(
      "Existing supplemental cache was generated with different settings; ",
      "recomputing fast panels and resuming any matching per-split caches."
    )
  }

  result_dir <- .get_results_dir()
  deseq_results <- .sup_read_result(
    deseq_results,
    file.path(result_dir, "analysis", "deseq", "deseq_results.C1.rds"),
    "DESeq2"
  )
  gene_se <- .sup_read_result(
    gene_se,
    file.path(PROJECT_ROOT, "resources", "cache", "tximeta_gse.rds"),
    "gene-level tximeta SE",
    required = "DESeq2" %in% permutation_tools && isTRUE(run_permutations)
  )
  dexseq_results <- .sup_read_result(
    dexseq_results,
    .sup_find_dexseq_cache(result_dir),
    "DEXSeq",
    required = "DEXSeq" %in% permutation_tools && isTRUE(run_permutations)
  )
  drimseq_results <- .sup_read_result(
    drimseq_results,
    file.path(result_dir, "analysis", "drimseq", "drimseq_results.C1.rds"),
    "DRIMSeq",
    required = "DRIMSeq" %in% permutation_tools && isTRUE(run_permutations)
  )
  suppa_results <- .sup_read_result(
    suppa_results,
    file.path(
      result_dir,
      "analysis", "suppa", "standard", "C1", "suppa_results.rds"
    ),
    "SUPPA2",
    required = "SUPPA2" %in% permutation_tools && isTRUE(run_permutations)
  )
  transcript_se <- .sup_read_result(
    transcript_se,
    file.path(PROJECT_ROOT, "resources", "cache", "tximeta_se.rds"),
    "transcript-level tximeta SE"
  )
  tx2gene <- .sup_read_result(
    tx2gene,
    file.path(PROJECT_ROOT, "resources", "cache", "tx2gene.rds"),
    "tx2gene"
  )

  message("Running gene-level PCA feature-set sensitivity analysis...")
  gene_pca_analysis <- .sup_run_gene_pca(
    deseq_results,
    feature_numbers = pca_feature_numbers
  )
  if (!(gene_pca_display %in% names(gene_pca_analysis$results))) {
    stop(
      "gene_pca_display not found: ", gene_pca_display,
      "\nAvailable: ", paste(names(gene_pca_analysis$results), collapse = ", ")
    )
  }
  gene_pca <- gene_pca_analysis$results[[gene_pca_display]]

  message("Running transcript-level PCA feature-set sensitivity analysis...")
  transcript_pca_analysis <- .sup_run_transcript_pca(
    transcript_se = transcript_se,
    tx2gene = tx2gene,
    tested_gene_ids = gene_pca_analysis$tested_gene_ids,
    feature_numbers = pca_feature_numbers,
    min_tpm = transcript_min_tpm,
    min_samples = transcript_min_samples,
    pseudocount = transcript_pseudocount
  )
  if (!(transcript_pca_display %in% names(transcript_pca_analysis$results))) {
    stop(
      "transcript_pca_display not found: ", transcript_pca_display,
      "\nAvailable: ",
      paste(names(transcript_pca_analysis$results), collapse = ", ")
    )
  }
  transcript_pca <- transcript_pca_analysis$results[[transcript_pca_display]]

  message("Running HALLMARK_HYPOXIA GSEA...")
  gsea <- .sup_run_hypoxia_gsea(
    deseq_results = deseq_results,
    transcript_se = transcript_se,
    tx2gene = tx2gene,
    hypoxia_genes = hypoxia_genes,
    include_transcript_ranking = include_transcript_gsea
  )

  pairwise_psi <- NULL
  if (isTRUE(run_pairwise_psi)) {
    message("Summarizing pairwise sample-level SUPPA2 PSI differences...")
    pairwise_psi <- .sup_run_pairwise_psi(
      deseq_results = deseq_results,
      psi_dir = suppa_psi_dir,
      cutoff = pairwise_psi_cutoff
    )
  }

  permutations <- NULL
  if (isTRUE(run_permutations)) {
    message(
      "Running C1/H1 balanced-label permutations. ",
      "Cached count matrices will be reused; label-dependent models will be refitted. ",
      "DRIMSeq precision mode: ", drimseq_precision_mode, "."
    )
    permutations <- .sup_run_permutations(
      deseq_results = deseq_results,
      gene_se = gene_se,
      dexseq_results = dexseq_results,
      drimseq_results = drimseq_results,
      suppa_results = suppa_results,
      permutation_tools = permutation_tools,
      n_null_splits = n_null_splits,
      alpha = padj_cutoff,
      suppa_threshold = suppa_threshold,
      suppa_master_tpm_file = suppa_master_tpm_file,
      suppa_ioe_file = suppa_ioe_file,
      suppa_executable = suppa_executable,
      suppa_runner = suppa_runner,
      drimseq_precision_mode = drimseq_precision_mode,
      workers = workers,
      seed = seed,
      cache_dir = file.path(outdir, "cache", "permutations"),
      force = force
    )
  }

  out <- list(
    meta = list(
      analysis_version = .SUP_FIG_X_ANALYSIS_VERSION,
      outdir = outdir,
      pca_feature_numbers = pca_feature_numbers,
      gene_pca_display = gene_pca_display,
      transcript_pca_display = transcript_pca_display,
      transcript_min_tpm = transcript_min_tpm,
      transcript_min_samples = transcript_min_samples,
      transcript_pseudocount = transcript_pseudocount,
      include_transcript_gsea = include_transcript_gsea,
      run_pairwise_psi = run_pairwise_psi,
      pairwise_psi_cutoff = pairwise_psi_cutoff,
      suppa_psi_dir = suppa_psi_dir,
      run_permutations = run_permutations,
      permutation_tools = permutation_tools,
      n_null_splits = n_null_splits,
      suppa_threshold = suppa_threshold,
      suppa_executable = suppa_executable,
      suppa_runner = suppa_runner,
      drimseq_precision_mode = drimseq_precision_mode,
      padj_cutoff = padj_cutoff,
      seed = seed
    ),
    pca = list(
      gene = gene_pca,
      transcript = transcript_pca,
      gene_sensitivity = gene_pca_analysis$results,
      transcript_sensitivity = transcript_pca_analysis$results,
      gene_meta = gene_pca_analysis[
        setdiff(names(gene_pca_analysis), "results")
      ],
      transcript_meta = transcript_pca_analysis[
        setdiff(names(transcript_pca_analysis), "results")
      ]
    ),
    gsea = gsea,
    pairwise_psi = pairwise_psi,
    permutations = permutations
  )

  if (isTRUE(save_tables)) {
    .sup_write_tsv(
      gene_pca$scores,
      file.path(outdir, "sup_fig_x_A_gene_pca_scores.tsv")
    )
    .sup_write_tsv(
      transcript_pca$scores,
      file.path(outdir, "sup_fig_x_B_transcript_pca_scores.tsv")
    )

    pca_sensitivity_scores <- do.call(rbind, c(
      lapply(names(gene_pca_analysis$results), function(feature_set) {
        transform(
          gene_pca_analysis$results[[feature_set]]$scores,
          modality = "Gene",
          feature_set = feature_set
        )
      }),
      lapply(names(transcript_pca_analysis$results), function(feature_set) {
        transform(
          transcript_pca_analysis$results[[feature_set]]$scores,
          modality = "Transcript",
          feature_set = feature_set
        )
      })
    ))
    rownames(pca_sensitivity_scores) <- NULL
    .sup_write_tsv(
      pca_sensitivity_scores,
      file.path(outdir, "sup_fig_x_AB_pca_sensitivity_scores.tsv")
    )

    pca_sensitivity_summary <- do.call(rbind, c(
      lapply(names(gene_pca_analysis$results), function(feature_set) {
        z <- gene_pca_analysis$results[[feature_set]]
        data.frame(
          modality = "Gene",
          feature_set = feature_set,
          n_features = z$n_features,
          pc1_variance_percent = z$variance_explained[1],
          pc2_variance_percent = z$variance_explained[2],
          transformation = gene_pca_analysis$transformation,
          universe_definition = gene_pca_analysis$universe_definition,
          stringsAsFactors = FALSE
        )
      }),
      lapply(names(transcript_pca_analysis$results), function(feature_set) {
        z <- transcript_pca_analysis$results[[feature_set]]
        data.frame(
          modality = "Transcript",
          feature_set = feature_set,
          n_features = z$n_features,
          pc1_variance_percent = z$variance_explained[1],
          pc2_variance_percent = z$variance_explained[2],
          transformation = transcript_pca_analysis$transformation,
          universe_definition = transcript_pca_analysis$universe_definition,
          stringsAsFactors = FALSE
        )
      })
    ))
    rownames(pca_sensitivity_summary) <- NULL
    .sup_write_tsv(
      pca_sensitivity_summary,
      file.path(outdir, "sup_fig_x_AB_pca_sensitivity_summary.tsv")
    )

    gsea_summary <- gsea$summary
    if ("leadingEdge" %in% colnames(gsea_summary)) {
      gsea_summary$leadingEdge <- vapply(
        gsea_summary$leadingEdge,
        paste,
        collapse = ",",
        FUN.VALUE = character(1)
      )
    }
    .sup_write_tsv(
      gsea_summary,
      file.path(outdir, "sup_fig_x_C_hypoxia_gsea.tsv")
    )

    if (!is.null(pairwise_psi)) {
      .sup_write_tsv(
        pairwise_psi$pair_summary,
        file.path(outdir, "sup_fig_x_E_pairwise_PSI_summary.tsv")
      )
      .sup_write_tsv(
        pairwise_psi$heatmap,
        file.path(outdir, "sup_fig_x_F_pairwise_PSI_heatmap.tsv")
      )
      .sup_write_tsv(
        pairwise_psi$partition_summary,
        file.path(outdir, "sup_fig_x_E_PSI_label_partitions.tsv")
      )
      .sup_write_tsv(
        pairwise_psi$partition_test,
        file.path(outdir, "sup_fig_x_E_PSI_partition_test.tsv")
      )
    }

    if (!is.null(permutations)) {
      .sup_write_tsv(
        permutations$membership,
        file.path(outdir, "sup_fig_x_D_permutation_membership.tsv")
      )
      .sup_write_tsv(
        permutations$tool_summary,
        file.path(outdir, "sup_fig_x_D_permutation_tool_summary.tsv")
      )
      .sup_write_tsv(
        permutations$plot_metrics,
        file.path(outdir, "sup_fig_x_D_permutation_plot_metrics.tsv")
      )
      .sup_write_tsv(
        permutations$split_descriptors,
        file.path(outdir, "sup_fig_x_D_permutation_split_key.tsv")
      )
      .sup_write_tsv(
        permutations$randomization_summary,
        file.path(outdir, "sup_fig_x_D_randomization_summary.tsv")
      )
      if (!is.null(permutations$deseq_observed_verification)) {
        .sup_write_tsv(
          permutations$deseq_observed_verification,
          file.path(outdir, "sup_fig_x_D_deseq_observed_verification.tsv")
        )
      }
    }
  }

  saveRDS(out, cache_file)
  out
}

# --------------------------------------------------
# Plotting functions
# --------------------------------------------------

.sup_pca_envelope_data <- function(
  df,
  expand = 1.18,
  min_radius = 0.035,
  n_points = 180L
) {
  if (
    length(expand) != 1 ||
      !is.finite(expand) ||
      expand <= 1
  ) {
    stop("PCA ellipse expansion must be one finite number greater than 1.")
  }
  if (
    length(min_radius) != 1 ||
      !is.finite(min_radius) ||
      min_radius <= 0
  ) {
    stop("PCA ellipse minimum radius must be a positive finite number.")
  }

  x_limits <- range(df$PC1, na.rm = TRUE)
  y_limits <- range(df$PC2, na.rm = TRUE)
  x_span <- diff(x_limits)
  y_span <- diff(y_limits)
  if (!is.finite(x_span) || x_span <= 0) x_span <- 1
  if (!is.finite(y_span) || y_span <= 0) y_span <- 1
  x_mid <- mean(x_limits)
  y_mid <- mean(y_limits)

  normalized <- data.frame(
    x = (df$PC1 - x_mid) / x_span,
    y = (df$PC2 - y_mid) / y_span,
    condition = df$condition,
    stringsAsFactors = FALSE
  )
  groups <- split(normalized, normalized$condition, drop = TRUE)
  groups <- groups[vapply(groups, nrow, integer(1)) >= 2]
  theta <- seq(0, 2 * pi, length.out = as.integer(n_points))

  envelopes <- lapply(names(groups), function(condition_name) {
    group_df <- groups[[condition_name]]
    coordinates <- as.matrix(group_df[, c("x", "y"), drop = FALSE])
    center <- colMeans(coordinates)
    covariance <- stats::cov(coordinates)
    if (
      any(!is.finite(covariance)) ||
        all(abs(covariance) < .Machine$double.eps)
    ) {
      covariance <- diag(min_radius^2, 2)
    }
    eig <- eigen(covariance, symmetric = TRUE)
    base_radii <- pmax(
      sqrt(pmax(eig$values, 0)),
      min_radius / expand
    )
    projected <- sweep(coordinates, 2, center, "-") %*% eig$vectors
    standardized_distance <- sqrt(rowSums(
      sweep(projected, 2, base_radii, "/")^2
    ))
    scale_factor <- max(standardized_distance, na.rm = TRUE) * expand
    if (!is.finite(scale_factor) || scale_factor <= 0) {
      scale_factor <- expand
    }
    radii <- pmax(base_radii * scale_factor, min_radius)
    ellipse_normalized <- cbind(
      radii[[1]] * cos(theta),
      radii[[2]] * sin(theta)
    ) %*% t(eig$vectors)
    ellipse_normalized <- sweep(
      ellipse_normalized,
      2,
      center,
      "+"
    )

    data.frame(
      PC1 = ellipse_normalized[, 1] * x_span + x_mid,
      PC2 = ellipse_normalized[, 2] * y_span + y_mid,
      condition = condition_name,
      stringsAsFactors = FALSE
    )
  })

  if (length(envelopes) == 0) {
    return(data.frame(
      PC1 = numeric(),
      PC2 = numeric(),
      condition = character()
    ))
  }
  do.call(rbind, envelopes)
}

.sup_plot_pca <- function(
  pca_result,
  title,
  base_size = 13,
  label_samples = TRUE,
  show_ellipses = FALSE,
  ellipse_expand = 1.18
) {
  .sup_check_pkg("ggplot2")
  df <- pca_result$scores
  variance <- pca_result$variance_explained
  colors <- .sup_condition_colors()
  df$condition <- factor(
    df$condition,
    levels = intersect(names(.sup_condition_labels), unique(df$condition))
  )

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = PC1, y = PC2, fill = condition)
  )

  if (isTRUE(show_ellipses)) {
    ellipse_df <- .sup_pca_envelope_data(
      df,
      expand = ellipse_expand
    )
    if (nrow(ellipse_df) > 0) {
      p <- p +
        ggplot2::geom_polygon(
          data = ellipse_df,
          ggplot2::aes(
            x = PC1,
            y = PC2,
            group = condition,
            color = condition,
            fill = condition
          ),
          inherit.aes = FALSE,
          alpha = 0.08,
          linewidth = 0.65,
          show.legend = FALSE
        ) +
        ggplot2::scale_color_manual(
          values = colors,
          breaks = names(.sup_condition_labels),
          drop = FALSE,
          guide = "none"
        )
    }
  }

  p <- p +
    ggplot2::geom_point(
      shape = 21,
      size = 4,
      stroke = 0.6,
      color = "black"
    ) +
    ggplot2::scale_fill_manual(
      values = colors,
      breaks = names(.sup_condition_labels),
      labels = .sup_condition_labels,
      drop = FALSE,
      name = NULL
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste0(
        pca_result$selection_label,
        " (n = ", format(pca_result$n_features, big.mark = ","), ")"
      ),
      x = sprintf("PC1 (%.1f%%)", variance[1]),
      y = sprintf("PC2 (%.1f%%)", variance[2])
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "vertical"
    )

  if (isTRUE(label_samples)) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        ggplot2::aes(label = sample_id),
        size = 3,
        box.padding = 0.35,
        point.padding = 0.35,
        max.overlaps = Inf,
        segment.color = NA,
        show.legend = FALSE,
        seed = 333
      )
    } else {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = sample_id),
        size = 3,
        vjust = -0.9,
        show.legend = FALSE
      )
    }
  }
  p
}

plot_sup_fig_x_pca_sensitivity <- function(
  res,
  outdir = NULL,
  base_size = 11
) {
  .sup_check_pkg(c("ggplot2", "patchwork"))
  outdir <- outdir %||% res$meta$outdir

  gene_plots <- lapply(res$pca$gene_sensitivity, function(z) {
    .sup_plot_pca(
      z,
      title = paste0("Gene: ", z$selection_label),
      base_size = base_size,
      label_samples = FALSE
    ) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(subtitle = paste0("n = ", format(z$n_features, big.mark = ",")))
  })

  transcript_plots <- lapply(res$pca$transcript_sensitivity, function(z) {
    .sup_plot_pca(
      z,
      title = paste0("Transcript: ", z$selection_label),
      base_size = base_size,
      label_samples = FALSE
    ) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(subtitle = paste0("n = ", format(z$n_features, big.mark = ",")))
  })

  combined <- patchwork::wrap_plots(
    c(gene_plots, transcript_plots),
    ncol = length(gene_plots),
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = "PCA Feature-Set Sensitivity",
      subtitle = paste(
        "Top-variable feature sets are compared with the complete",
        "DESeq2-tested gene and expressed-transcript universes."
      )
    ) &
    ggplot2::theme(legend.position = "bottom")

  ggplot2::ggsave(
    file.path(outdir, "sup_fig_x_AB_pca_feature_sensitivity.pdf"),
    combined,
    width = 16,
    height = 9
  )
  combined
}

plot_sup_fig_x_gene_pca <- function(
  res,
  outdir = NULL,
  base_size = 13,
  show_ellipses = FALSE,
  ellipse_expand = 1.18
) {
  outdir <- outdir %||% res$meta$outdir
  p <- .sup_plot_pca(
    res$pca$gene,
    title = "Gene-Level Sample PCA",
    base_size = base_size,
    show_ellipses = show_ellipses,
    ellipse_expand = ellipse_expand
  )
  filename <- if (isTRUE(show_ellipses)) {
    "sup_fig_x_A_gene_pca_ellipses.pdf"
  } else {
    "sup_fig_x_A_gene_pca.pdf"
  }
  ggplot2::ggsave(
    file.path(outdir, filename),
    p,
    width = 7,
    height = 5.5
  )
  p
}

plot_sup_fig_x_transcript_pca <- function(
  res,
  outdir = NULL,
  base_size = 13,
  show_ellipses = FALSE,
  ellipse_expand = 1.18
) {
  outdir <- outdir %||% res$meta$outdir
  p <- .sup_plot_pca(
    res$pca$transcript,
    title = "Transcript-Level Sample PCA",
    base_size = base_size,
    show_ellipses = show_ellipses,
    ellipse_expand = ellipse_expand
  )
  filename <- if (isTRUE(show_ellipses)) {
    "sup_fig_x_B_transcript_pca_ellipses.pdf"
  } else {
    "sup_fig_x_B_transcript_pca.pdf"
  }
  ggplot2::ggsave(
    file.path(outdir, filename),
    p,
    width = 7,
    height = 5.5
  )
  p
}

plot_sup_fig_x_hypoxia_gsea <- function(
  res,
  outdir = NULL,
  base_size = 13
) {
  .sup_check_pkg("ggplot2")
  outdir <- outdir %||% res$meta$outdir

  gsea_parts <- res$gsea[
    intersect(c("gene", "transcript_derived"), names(res$gsea))
  ]
  running <- do.call(rbind, lapply(gsea_parts, `[[`, "running"))
  summary <- res$gsea$summary

  annotation <- data.frame(
    ranking = summary$ranking,
    label = sprintf(
      "NES = %.2f\nFDR = %.3g",
      summary$NES,
      summary$padj
    ),
    stringsAsFactors = FALSE
  )
  annotation$x <- vapply(split(running$rank, running$ranking), max, numeric(1))[
    annotation$ranking
  ] * 0.98
  annotation$y <- vapply(
    split(running$running_score, running$ranking),
    max,
    numeric(1),
    na.rm = TRUE
  )[annotation$ranking]

  p <- ggplot2::ggplot(
    running,
    ggplot2::aes(x = rank, y = running_score)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey70", linewidth = 0.4) +
    ggplot2::geom_line(color = "#0072B2", linewidth = 0.8) +
    ggplot2::geom_rug(
      data = running[running$hit, , drop = FALSE],
      sides = "b",
      alpha = 0.25,
      linewidth = 0.25
    ) +
    ggplot2::geom_label(
      data = annotation,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      size = 3,
      label.size = 0.2
    ) +
    ggplot2::facet_wrap(~ranking, ncol = 1, scales = "free_x") +
    ggplot2::labs(
      title = "Hallmark Hypoxia Gene-Set Enrichment",
      x = "Ranked Feature Position",
      y = "Running Enrichment Score"
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )

  ggplot2::ggsave(
    file.path(outdir, "sup_fig_x_C_hypoxia_gsea.pdf"),
    p,
    width = 7,
    height = if (length(gsea_parts) > 1) 7 else 4.5
  )
  p
}

plot_sup_fig_x_permutations <- function(
  res,
  outdir = NULL,
  base_size = 13
) {
  .sup_check_pkg(c("ggplot2", "scales", "patchwork"))
  if (is.null(res$permutations)) {
    stop("Permutation analysis was not run.")
  }
  outdir <- outdir %||% res$meta$outdir

  df <- res$permutations$plot_metrics
  keep_metrics <- c("DESeq2", "DEXSeq", "DRIMSeq", "SUPPA2")
  df <- df[df$metric %in% keep_metrics, , drop = FALSE]
  active_metrics <- keep_metrics[
    keep_metrics %in% unique(df$metric[!is.na(df$n_significant)])
  ]
  if (length(active_metrics) == 0) {
    stop("No completed permutation metrics are available to plot.")
  }
  df <- df[df$metric %in% active_metrics, , drop = FALSE]
  df$metric <- factor(df$metric, levels = active_metrics)

  split_key <- res$permutations$split_descriptors
  if (is.null(split_key)) {
    split_key <- .sup_split_descriptors(res$permutations$membership)
  }
  is_observed_split <- split_key$split_short == "Observed"
  split_key$split_number <- integer(nrow(split_key))
  split_key$split_number[is_observed_split] <- 0L
  split_key$split_number[!is_observed_split] <- as.integer(
    sub("^P", "", split_key$split_short[!is_observed_split])
  )
  split_levels <- split_key$split_short[order(split_key$split_number)]
  df$split_short <- factor(df$split_short, levels = split_levels)

  colors <- .sup_tool_colors()
  tool_plots <- lapply(active_metrics, function(metric_name) {
    tool_df <- df[df$metric == metric_name, , drop = FALSE]
    label_df <- tool_df[
      is.finite(tool_df$n_significant),
      ,
      drop = FALSE
    ]
    null_df <- tool_df[
      tool_df$split_type == "Permuted",
      ,
      drop = FALSE
    ]
    observed_df <- tool_df[
      tool_df$split_type == "Observed",
      ,
      drop = FALSE
    ]
    tool_color <- unname(colors[[metric_name]])

    ggplot2::ggplot(
      tool_df,
      ggplot2::aes(x = split_short, y = n_significant)
    ) +
      ggplot2::geom_point(
        data = null_df,
        shape = 21,
        size = 3.2,
        stroke = 0.45,
        fill = tool_color,
        color = "black",
        alpha = 0.72
      ) +
      ggplot2::geom_point(
        data = observed_df,
        shape = 23,
        size = 4.8,
        stroke = 0.7,
        fill = tool_color,
        color = "black"
      ) +
      ggplot2::geom_text(
        data = label_df,
        ggplot2::aes(label = scales::comma(n_significant)),
        vjust = -1.0,
        size = 2.8,
        show.legend = FALSE,
        na.rm = TRUE
      ) +
      ggplot2::scale_x_discrete(drop = FALSE) +
      ggplot2::scale_y_continuous(
        labels = scales::label_comma(),
        expand = ggplot2::expansion(mult = c(0.04, 0.26))
      ) +
      ggplot2::labs(
        title = metric_name,
        x = "Label Assignment",
        y = "Significant Genes"
      ) +
      ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
        plot.title = ggplot2::element_text(face = "bold"),
        plot.margin = ggplot2::margin(t = 8, r = 8, b = 8, l = 8)
      )
  })
  names(tool_plots) <- active_metrics

  permutation_subtitle <- paste0(
    "Diamonds show the observed labels; circles show ",
    res$permutations$meta$n_null_splits_run,
    if (isTRUE(res$permutations$meta$exhaustive)) {
      " exhaustive alternative balanced assignments"
    } else {
      " pilot alternative balanced assignments"
    },
    ".\nValues label every assignment",
    if ("SUPPA2" %in% active_metrics) {
      paste0(
        ". SUPPA2 uses ",
        if (res$permutations$meta$suppa_threshold == "paper_raw_p") {
          "the manuscript threshold (unadjusted P)"
        } else {
          "BH-adjusted P values"
        }
      )
    } else {
      ""
    },
    "."
  )

  p_grid <- patchwork::wrap_plots(
    tool_plots,
    ncol = 2
  ) +
    patchwork::plot_annotation(
      title = "D. C1/H1 Balanced-Label Diagnostic",
      subtitle = permutation_subtitle,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        plot.subtitle = ggplot2::element_text(
          size = ggplot2::rel(0.88)
        )
      )
    )

  key_df <- split_key[!is_observed_split, , drop = FALSE]
  key_df <- key_df[order(key_df$split_number), , drop = FALSE]
  key_df$key_label <- paste0(
    key_df$split_short,
    ": ",
    key_df$swapped_C1_sample,
    "  \u2194  ",
    key_df$swapped_H1_sample
  )
  key_ncol <- min(3L, nrow(key_df))
  key_nrow <- ceiling(nrow(key_df) / key_ncol)
  key_index <- seq_len(nrow(key_df))
  key_df$key_col <- ((key_index - 1L) %/% key_nrow) + 1L
  key_df$key_y <- key_nrow - ((key_index - 1L) %% key_nrow)

  p_key <- ggplot2::ggplot(
    key_df,
    ggplot2::aes(x = key_col, y = key_y, label = key_label)
  ) +
    ggplot2::geom_text(
      hjust = 0,
      size = max(2.8, base_size / 4.5),
      lineheight = 1.05
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0.75, key_ncol + 0.8),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0.4, key_nrow + 0.7),
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = "Label-Assignment Key",
      subtitle = paste(
        "P labels denote the same C1-H1 sample swap in every tool;",
        "diamond = observed labels, circle = alternative assignment."
      )
    ) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(
        size = ggplot2::rel(0.85),
        margin = ggplot2::margin(b = 8)
      ),
      plot.margin = ggplot2::margin(t = 8, r = 18, b = 8, l = 18)
    )

  p_saved <- patchwork::wrap_plots(
    p_grid,
    p_key,
    ncol = 1,
    heights = c(4.8, 1.55)
  )

  ggplot2::ggsave(
    file.path(outdir, "sup_fig_x_D_label_permutations_v2.pdf"),
    p_saved,
    width = 14,
    height = if (length(active_metrics) > 2) 11.5 else 8.2
  )
  p_saved
}

plot_sup_fig_x_pairwise_psi <- function(
  res,
  outdir = NULL,
  base_size = 13
) {
  .sup_check_pkg(c("ggplot2", "scales", "patchwork"))
  if (is.null(res$pairwise_psi)) {
    stop("Pairwise SUPPA2 PSI analysis was not run.")
  }
  outdir <- outdir %||% res$meta$outdir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  partition_df <- res$pairwise_psi$partition_summary
  partition_test <- res$pairwise_psi$partition_test
  split_number <- integer(nrow(partition_df))
  is_observed <- partition_df$split_type == "Observed"
  split_number[is_observed] <- 0L
  split_number[!is_observed] <- as.integer(
    sub("^P", "", partition_df$split_short[!is_observed])
  )
  split_levels <- partition_df$split_short[order(split_number)]
  partition_df$split_short <- factor(
    partition_df$split_short,
    levels = split_levels
  )
  observed_df <- partition_df[is_observed, , drop = FALSE]
  null_df <- partition_df[!is_observed, , drop = FALSE]

  p_partition <- ggplot2::ggplot(
    partition_df,
    ggplot2::aes(
      x = split_short,
      y = between_minus_within_percent_points
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "grey50",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggplot2::geom_point(
      data = null_df,
      shape = 21,
      size = 3.2,
      fill = "#7A5195",
      color = "#4A2E61",
      stroke = 0.5,
      alpha = 0.8
    ) +
    ggplot2::geom_point(
      data = observed_df,
      shape = 23,
      size = 4.8,
      fill = "#B2182B",
      color = "black",
      stroke = 0.7
    ) +
    ggplot2::geom_text(
      data = observed_df,
      ggplot2::aes(
        label = paste0(
          scales::number(
            between_minus_within_percent_points,
            accuracy = 0.01
          ),
          " pp"
        )
      ),
      vjust = -1.1,
      size = 3.2
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_number(accuracy = 0.1),
      expand = ggplot2::expansion(mult = c(0.08, 0.18))
    ) +
    ggplot2::labs(
      title = "Exact Label-Permutation Test of PSI Separation",
      subtitle = paste0(
        "Diamond: observed labels; circles: all nine alternatives. ",
        "Observed rank ",
        partition_test$observed_separation_rank,
        "/",
        partition_test$n_unique_balanced_partitions,
        "; exact P = ",
        scales::number(partition_test$exact_upper_tail_p, accuracy = 0.01)
      ),
      caption = paste0(
        "Distance = percentage of ",
        scales::comma(res$pairwise_psi$n_complete_events),
        " complete-case events with |Delta PSI| >= ",
        format(res$pairwise_psi$cutoff, nsmall = 2),
        ". Statistic = mean between-group distance - mean within-group distance."
      ),
      x = "Label Assignment",
      y = "Between-Group Minus Within-Group Distance\n(Percentage Points)"
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
      plot.caption = ggplot2::element_text(
        hjust = 0,
        size = ggplot2::rel(0.78)
      )
    )

  heat_df <- res$pairwise_psi$heatmap
  sample_order <- res$pairwise_psi$sample_order
  sample_display <- paste0(
    sample_order,
    "\n",
    ifelse(
      res$pairwise_psi$sample_condition[sample_order] == "C1",
      "Normoxia",
      "Hypoxia (1 h)"
    )
  )
  names(sample_display) <- sample_order
  heat_df$sample_1_index <- match(heat_df$sample_1, sample_order)
  heat_df$sample_2_index <- match(heat_df$sample_2, sample_order)
  heat_df <- heat_df[
    heat_df$sample_1_index <= heat_df$sample_2_index,
    ,
    drop = FALSE
  ]
  heat_df$is_diagonal <- heat_df$sample_1_index == heat_df$sample_2_index
  heat_df$plot_value <- heat_df$percent_events_abs_dpsi_ge_cutoff
  heat_df$plot_value[heat_df$is_diagonal] <- NA_real_
  heat_labels <- heat_df[!heat_df$is_diagonal, , drop = FALSE]
  fill_limits <- range(heat_labels$plot_value, na.rm = TRUE)
  if (!all(is.finite(fill_limits))) {
    stop("No finite off-diagonal PSI distances are available to plot.")
  }
  if (diff(fill_limits) == 0) {
    fill_limits <- fill_limits + c(-0.5, 0.5)
  }
  heat_df$sample_1 <- factor(heat_df$sample_1, levels = sample_order)
  heat_df$sample_2 <- factor(
    heat_df$sample_2,
    levels = rev(sample_order)
  )
  heat_labels$sample_1 <- factor(
    heat_labels$sample_1,
    levels = sample_order
  )
  heat_labels$sample_2 <- factor(
    heat_labels$sample_2,
    levels = rev(sample_order)
  )
  p_heat <- ggplot2::ggplot(
    heat_df,
    ggplot2::aes(
      x = sample_1,
      y = sample_2,
      fill = plot_value
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      data = heat_labels,
      ggplot2::aes(
        label = scales::number(
          plot_value,
          accuracy = 0.1
        )
      ),
      size = 3
    ) +
    ggplot2::scale_fill_viridis_c(
      option = "C",
      direction = -1,
      limits = fill_limits,
      na.value = "#E5E5E5",
      name = "% Events"
    ) +
    ggplot2::scale_x_discrete(labels = sample_display) +
    ggplot2::scale_y_discrete(labels = sample_display) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = "Pairwise PSI Difference Matrix",
      subtitle = paste0(
        scales::comma(res$pairwise_psi$n_complete_events),
        " events with finite PSI in all six samples; |Delta PSI| >= ",
        format(res$pairwise_psi$cutoff, nsmall = 2),
        ". Diagonal shown in grey."
      ),
      x = "Sample",
      y = "Sample"
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical"
    )

  p_combined <- patchwork::wrap_plots(
    p_partition + ggplot2::labs(tag = "E"),
    p_heat + ggplot2::labs(tag = "F"),
    ncol = 2,
    widths = c(1.08, 1)
  )

  ggplot2::ggsave(
    file.path(outdir, "sup_fig_x_E_pairwise_PSI_v2.pdf"),
    p_partition,
    width = 8.5,
    height = 5.8
  )
  ggplot2::ggsave(
    file.path(outdir, "sup_fig_x_F_pairwise_PSI_heatmap_v2.pdf"),
    p_heat,
    width = 7.5,
    height = 6.2
  )
  ggplot2::ggsave(
    file.path(outdir, "sup_fig_x_EF_pairwise_PSI_v2.pdf"),
    p_combined,
    width = 15,
    height = 6.5
  )

  list(
    dot = p_partition,
    partition = p_partition,
    heatmap = p_heat,
    combined = p_combined
  )
}

# --------------------------------------------------
# Convenience wrapper
# --------------------------------------------------

plot_sup_fig_x_all <- function(
  deseq_results = NULL,
  gene_se = NULL,
  dexseq_results = NULL,
  drimseq_results = NULL,
  suppa_results = NULL,
  transcript_se = NULL,
  tx2gene = NULL,
  outdir = NULL,
  pca_feature_numbers = c(500L, 1000L, 5000L),
  gene_pca_display = "all_tested",
  transcript_pca_display = "all_expressed",
  pca_ellipse_expand = 1.18,
  transcript_min_tpm = 1,
  transcript_min_samples = 3L,
  transcript_pseudocount = 0.5,
  include_transcript_gsea = TRUE,
  hypoxia_genes = NULL,
  run_pairwise_psi = TRUE,
  pairwise_psi_cutoff = 0.10,
  suppa_psi_dir = NULL,
  run_permutations = TRUE,
  permutation_tools = c("DESeq2", "DEXSeq", "DRIMSeq", "SUPPA2"),
  n_null_splits = Inf,
  padj_cutoff = 0.10,
  suppa_threshold = c("paper_raw_p", "BH"),
  suppa_master_tpm_file = NULL,
  suppa_ioe_file = NULL,
  suppa_executable = "suppa.py",
  suppa_runner = file.path(PROJECT_ROOT, "scripts", "run"),
  drimseq_precision_mode = c("reuse_observed", "refit"),
  workers = 1L,
  seed = 333L,
  save_tables = TRUE,
  force = FALSE,
  base_size = 13
) {
  .sup_check_pkg(c("ggplot2", "patchwork"))
  suppa_threshold <- match.arg(suppa_threshold)
  drimseq_precision_mode <- match.arg(drimseq_precision_mode)
  if (is.null(outdir)) {
    outdir <- file.path(.get_results_dir(), "plots", "sup_fig_x")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  res <- run_sup_fig_x_analysis(
    deseq_results = deseq_results,
    gene_se = gene_se,
    dexseq_results = dexseq_results,
    drimseq_results = drimseq_results,
    suppa_results = suppa_results,
    transcript_se = transcript_se,
    tx2gene = tx2gene,
    outdir = outdir,
    pca_feature_numbers = pca_feature_numbers,
    gene_pca_display = gene_pca_display,
    transcript_pca_display = transcript_pca_display,
    transcript_min_tpm = transcript_min_tpm,
    transcript_min_samples = transcript_min_samples,
    transcript_pseudocount = transcript_pseudocount,
    include_transcript_gsea = include_transcript_gsea,
    hypoxia_genes = hypoxia_genes,
    run_pairwise_psi = run_pairwise_psi,
    pairwise_psi_cutoff = pairwise_psi_cutoff,
    suppa_psi_dir = suppa_psi_dir,
    run_permutations = run_permutations,
    permutation_tools = permutation_tools,
    n_null_splits = n_null_splits,
    padj_cutoff = padj_cutoff,
    suppa_threshold = suppa_threshold,
    suppa_master_tpm_file = suppa_master_tpm_file,
    suppa_ioe_file = suppa_ioe_file,
    suppa_executable = suppa_executable,
    suppa_runner = suppa_runner,
    drimseq_precision_mode = drimseq_precision_mode,
    workers = workers,
    seed = seed,
    save_tables = save_tables,
    force = force
  )

  p_a <- plot_sup_fig_x_gene_pca(
    res,
    outdir,
    base_size,
    show_ellipses = FALSE,
    ellipse_expand = pca_ellipse_expand
  )
  p_b <- plot_sup_fig_x_transcript_pca(
    res,
    outdir,
    base_size,
    show_ellipses = FALSE,
    ellipse_expand = pca_ellipse_expand
  )
  p_a_ellipses <- plot_sup_fig_x_gene_pca(
    res,
    outdir,
    base_size,
    show_ellipses = TRUE,
    ellipse_expand = pca_ellipse_expand
  )
  p_b_ellipses <- plot_sup_fig_x_transcript_pca(
    res,
    outdir,
    base_size,
    show_ellipses = TRUE,
    ellipse_expand = pca_ellipse_expand
  )
  p_ab_sensitivity <- plot_sup_fig_x_pca_sensitivity(
    res,
    outdir,
    base_size = max(9, base_size - 2)
  )
  p_c <- plot_sup_fig_x_hypoxia_gsea(res, outdir, base_size)
  p_d <- if (isTRUE(run_permutations)) {
    plot_sup_fig_x_permutations(res, outdir, base_size)
  } else {
    patchwork::plot_spacer()
  }
  psi_plots <- if (isTRUE(run_pairwise_psi)) {
    plot_sup_fig_x_pairwise_psi(res, outdir, base_size)
  } else {
    list(
      dot = patchwork::plot_spacer(),
      heatmap = patchwork::plot_spacer(),
      combined = NULL
    )
  }
  p_e <- psi_plots$dot
  p_f <- psi_plots$heatmap

  assemble_combined <- function(pca_a, pca_b) {
    pca_column <- patchwork::wrap_plots(
      pca_a + ggplot2::labs(tag = "A"),
      pca_b + ggplot2::labs(tag = "B"),
      ncol = 1,
      heights = c(1, 1),
      guides = "collect"
    ) &
      ggplot2::theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        plot.tag = ggplot2::element_text(face = "bold", size = 14)
      )

    gsea_column <- p_c +
      ggplot2::labs(tag = "C") +
      ggplot2::theme(
        plot.tag = ggplot2::element_text(face = "bold", size = 14)
      )

    top_section <- patchwork::wrap_plots(
      pca_column,
      gsea_column,
      ncol = 2,
      widths = c(1, 1.05)
    )

    psi_section <- patchwork::wrap_plots(
      p_e +
        ggplot2::labs(tag = "E") +
        ggplot2::theme(
          plot.tag = ggplot2::element_text(face = "bold", size = 14)
        ),
      p_f +
        ggplot2::labs(tag = "F") +
        ggplot2::theme(
          legend.position = "right",
          legend.direction = "vertical",
          legend.box = "vertical",
          plot.tag = ggplot2::element_text(face = "bold", size = 14)
        ),
      ncol = 2,
      widths = c(1.08, 1)
    )

    patchwork::wrap_plots(
      top_section,
      p_d,
      psi_section,
      ncol = 1,
      heights = c(1.65, 2.35, 1.25)
    )
  }

  combined <- assemble_combined(p_a, p_b)
  combined_ellipses <- assemble_combined(
    p_a_ellipses,
    p_b_ellipses
  )

  ggplot2::ggsave(
    file.path(outdir, "sup_fig_x_combined_v2.pdf"),
    combined,
    width = 14,
    height = 31
  )
  ggplot2::ggsave(
    file.path(outdir, "sup_fig_x_combined_pca_ellipses_v2.pdf"),
    combined_ellipses,
    width = 14,
    height = 31
  )

  invisible(list(
    analysis = res,
    fig_x_a = p_a,
    fig_x_b = p_b,
    fig_x_a_ellipses = p_a_ellipses,
    fig_x_b_ellipses = p_b_ellipses,
    pca_sensitivity = p_ab_sensitivity,
    fig_x_c = p_c,
    fig_x_d = p_d,
    fig_x_e = p_e,
    fig_x_f = p_f,
    pairwise_psi = psi_plots$combined,
    combined = combined,
    combined_pca_ellipses = combined_ellipses
  ))
}

# --------------------------------------------------
# Suggested calls
# --------------------------------------------------
#
# Full reviewer-facing analysis (default; all nine alternative splits):
#
#   source(file.path(PROJECT_ROOT, "scripts/plots/plot_sup_fig_x_v2.R"))
#   sup_fig_x <- plot_sup_fig_x_all(
#     drimseq_precision_mode = "reuse_observed",
#     workers = 2
#   )
#
# The wrapper writes both sup_fig_x_combined_v2.pdf (repelled sample labels
# without connector segments) and sup_fig_x_combined_pca_ellipses_v2.pdf (the
# same figure with descriptive condition envelopes in panels A and B). The
# envelope padding can be adjusted with pca_ellipse_expand; larger values
# produce wider ovals. The PCA legend is collected only within the upper-left
# PCA column; the PSI heatmap legend remains vertical and to the right.
#
# The pairwise PSI panels use the existing standard SUPPA2 files
# events/events_C1.psi and events/events_H1.psi. Completed per-tool
# permutation fits are resumed from outdir/cache/permutations when force=FALSE
# (the default), including fits created by the preceding script version.
# PSI distances use the same complete-case event universe for every sample pair.
# Panel E exhaustively evaluates all 10 unique balanced C1/H1 partitions; this
# fast sample-distance test is independent of n_null_splits.
# SUPPA2 permutation commands default to:
#   "$PROJECT_ROOT/scripts/run" suppa.py <subcommand> ...
# Set suppa_runner = NULL only when suppa.py is already available on the R PATH.
# On Windows, two workers are usually safer than four for DRIMSeq because
# SnowParam starts separate R processes and can multiply memory use.
#
# To display one variable-feature PCA in the combined figure instead:
#
#   sup_fig_x_top1000 <- plot_sup_fig_x_all(
#     gene_pca_display = "top_1000",
#     transcript_pca_display = "top_1000",
#     workers = 2
#   )
#
# Fast smoke test before committing compute:
#
#   sup_fig_x_pilot <- plot_sup_fig_x_all(
#     permutation_tools = "DESeq2",
#     n_null_splits = 1,
#     outdir = file.path(.get_results_dir(), "plots", "sup_fig_x_pilot"),
#     force = TRUE
#   )
#
# A one-shuffle result is intentionally labeled as a pilot. For the manuscript,
# run all nine alternative splits (n_null_splits = Inf).
