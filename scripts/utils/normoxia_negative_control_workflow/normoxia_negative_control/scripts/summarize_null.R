#!/usr/bin/env Rscript

parse_args <- function(x) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    if (!startsWith(x[[i]], "--") || i == length(x)) {
      stop("Arguments must be supplied as --name value pairs.")
    }
    out[[sub("^--", "", x[[i]])]] <- x[[i + 1L]]
    i <- i + 2L
  }
  out
}

`%||%` <- function(x, y) if (is.null(x)) y else x
args <- parse_args(commandArgs(trailingOnly = TRUE))
analysis_root <- normalizePath(
  args[["analysis-root"]] %||% stop("--analysis-root is required."),
  mustWork = TRUE
)
alpha <- as.numeric(args[["alpha"]] %||% 0.10)
if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) stop("--alpha must lie between 0 and 1.")

write_tsv <- function(x, path) {
  utils::write.table(
    x, file = path, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = ""
  )
}

strip_version <- function(x) sub("\\.[0-9]+$", "", as.character(x))

statistics_path <- file.path(analysis_root, "statistics_summary.tsv")
membership_path <- file.path(analysis_root, "permutation_membership.tsv")
sig_gene_path <- file.path(analysis_root, "significant_gene_ids.tsv")
if (!file.exists(statistics_path)) stop("Missing: ", statistics_path)
if (!file.exists(membership_path)) stop("Missing: ", membership_path)

statistics <- utils::read.delim(
  statistics_path, stringsAsFactors = FALSE, check.names = FALSE
)
membership <- utils::read.delim(
  membership_path, stringsAsFactors = FALSE, check.names = FALSE
)
if (file.exists(sig_gene_path)) {
  significant_genes <- utils::read.delim(
    sig_gene_path, stringsAsFactors = FALSE, check.names = FALSE
  )
} else {
  significant_genes <- data.frame(
    perm_id = character(), tool = character(), gene_id = character(),
    stringsAsFactors = FALSE
  )
}

suppa_rows <- list()
suppa_sig_rows <- list()
row_i <- 0L
sig_i <- 0L

add_suppa_row <- function(perm_id, tool, level, tested, significant, threshold, status = "ok") {
  row_i <<- row_i + 1L
  suppa_rows[[row_i]] <<- data.frame(
    perm_id = perm_id,
    tool = tool,
    level = level,
    n_tested = as.integer(tested),
    n_significant = as.integer(significant),
    fraction_significant = if (is.finite(tested) && tested > 0) significant / tested else NA_real_,
    threshold = threshold,
    status = status,
    stringsAsFactors = FALSE
  )
}

perm_ids <- unique(membership$perm_id)
for (perm_id in perm_ids) {
  dpsi_path <- file.path(analysis_root, "suppa", perm_id, "diffSplice.dpsi")
  if (!file.exists(dpsi_path)) {
    for (tool in c("SUPPA2_rawP", "SUPPA2_BH")) {
      add_suppa_row(perm_id, tool, "event", NA, NA, "missing", "dpsi file missing")
      add_suppa_row(perm_id, tool, "gene", NA, NA, "missing", "dpsi file missing")
    }
    next
  }

  x <- utils::read.delim(
    dpsi_path,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (ncol(x) < 2L) {
    stop("Expected at least two columns in: ", dpsi_path)
  }
  event_id <- rownames(x)
  cn <- colnames(x)
  p_col <- grep("p[-._]?val|pvalue", cn, ignore.case = TRUE, value = TRUE)
  if (!length(p_col)) p_col <- cn[[2]]
  pvalue <- suppressWarnings(as.numeric(x[[p_col[[1]]]]))
  padj <- stats::p.adjust(pvalue, method = "BH")
  gene_id <- strip_version(sub(";.*$", "", event_id))
  tested <- !is.na(pvalue)
  raw_sig <- tested & pvalue < alpha
  bh_sig <- tested & padj < alpha

  add_suppa_row(
    perm_id, "SUPPA2_rawP", "event",
    sum(tested), sum(raw_sig),
    sprintf("unadjusted p < %.3g", alpha)
  )
  add_suppa_row(
    perm_id, "SUPPA2_rawP", "gene",
    length(unique(gene_id[tested])), length(unique(gene_id[raw_sig])),
    sprintf("any event unadjusted p < %.3g", alpha)
  )
  add_suppa_row(
    perm_id, "SUPPA2_BH", "event",
    sum(tested), sum(bh_sig),
    sprintf("BH FDR < %.3g", alpha)
  )
  add_suppa_row(
    perm_id, "SUPPA2_BH", "gene",
    length(unique(gene_id[tested])), length(unique(gene_id[bh_sig])),
    sprintf("any event BH FDR < %.3g", alpha)
  )

  bh_genes <- sort(unique(gene_id[bh_sig]))
  if (length(bh_genes)) {
    sig_i <- sig_i + 1L
    suppa_sig_rows[[sig_i]] <- data.frame(
      perm_id = perm_id,
      tool = "SUPPA2_BH",
      gene_id = bh_genes,
      stringsAsFactors = FALSE
    )
  }

  event_results <- data.frame(
    event_id = event_id,
    gene_id = gene_id,
    pvalue = pvalue,
    padj_BH = padj,
    significant_rawP = raw_sig,
    significant_BH = bh_sig,
    stringsAsFactors = FALSE
  )
  write_tsv(
    event_results,
    file.path(analysis_root, "suppa", perm_id, "diffSplice_with_BH.tsv")
  )
}

suppa_summary <- do.call(rbind, suppa_rows)
combined <- rbind(statistics, suppa_summary)
write_tsv(combined, file.path(analysis_root, "all_tool_results_by_permutation.tsv"))

if (length(suppa_sig_rows)) {
  suppa_sig <- do.call(rbind, suppa_sig_rows)
  significant_genes <- rbind(significant_genes, suppa_sig)
}
significant_genes$gene_id <- strip_version(significant_genes$gene_id)
write_tsv(significant_genes, file.path(analysis_root, "all_significant_gene_ids.tsv"))

summarize_group <- function(z) {
  tool_name <- z$tool[[1]]
  level_name <- z$level[[1]]
  ok <- z$status == "ok" & !is.na(z$n_tested) & !is.na(z$n_significant)
  z <- z[ok, , drop = FALSE]
  if (!nrow(z)) {
    return(data.frame(
      tool = tool_name,
      level = level_name,
      n_permutations = 0L,
      median_tested = NA_real_,
      median_significant = NA_real_,
      q1_significant = NA_real_,
      q3_significant = NA_real_,
      min_significant = NA_real_,
      max_significant = NA_real_,
      median_fraction = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    tool = z$tool[[1]],
    level = z$level[[1]],
    n_permutations = nrow(z),
    median_tested = stats::median(z$n_tested),
    median_significant = stats::median(z$n_significant),
    q1_significant = unname(stats::quantile(z$n_significant, 0.25)),
    q3_significant = unname(stats::quantile(z$n_significant, 0.75)),
    min_significant = min(z$n_significant),
    max_significant = max(z$n_significant),
    median_fraction = stats::median(z$fraction_significant),
    stringsAsFactors = FALSE
  )
}

group_key <- interaction(combined$tool, combined$level, drop = TRUE, lex.order = TRUE)
aggregate_summary <- do.call(rbind, lapply(split(combined, group_key), summarize_group))
rownames(aggregate_summary) <- NULL
write_tsv(aggregate_summary, file.path(analysis_root, "all_tool_null_summary.tsv"))

get_count <- function(perm_id, tool) {
  z <- combined[
    combined$perm_id == perm_id &
      combined$tool == tool &
      combined$level == "gene" &
      combined$status == "ok",
    "n_significant"
  ]
  if (length(z) == 1L) as.integer(z) else NA_integer_
}

get_sig_ids <- function(perm_id, tool) {
  unique(significant_genes$gene_id[
    significant_genes$perm_id == perm_id & significant_genes$tool == tool
  ])
}

gene_counts <- do.call(rbind, lapply(perm_ids, function(perm_id) {
  dex <- get_sig_ids(perm_id, "DEXSeq")
  dri <- get_sig_ids(perm_id, "DRIMSeq")
  sup <- get_sig_ids(perm_id, "SUPPA2_BH")
  three_way <- Reduce(intersect, list(dex, dri, sup))
  splicing_union <- Reduce(union, list(dex, dri, sup))
  data.frame(
    perm_id = perm_id,
    DESeq2_DEG = get_count(perm_id, "DESeq2"),
    DEXSeq_AS_gene = get_count(perm_id, "DEXSeq"),
    DRIMSeq_AS_gene = get_count(perm_id, "DRIMSeq"),
    SUPPA2_BH_AS_gene = get_count(perm_id, "SUPPA2_BH"),
    AS_three_way_intersection = length(three_way),
    AS_three_tool_union = length(splicing_union),
    stringsAsFactors = FALSE
  )
}))
write_tsv(gene_counts, file.path(analysis_root, "null_gene_counts_by_permutation.tsv"))

count_columns <- setdiff(colnames(gene_counts), "perm_id")
gene_count_summary <- do.call(rbind, lapply(count_columns, function(nm) {
  x <- gene_counts[[nm]]
  x <- x[!is.na(x)]
  data.frame(
    metric = nm,
    n_permutations = length(x),
    median = if (length(x)) stats::median(x) else NA_real_,
    q1 = if (length(x)) unname(stats::quantile(x, 0.25)) else NA_real_,
    q3 = if (length(x)) unname(stats::quantile(x, 0.75)) else NA_real_,
    min = if (length(x)) min(x) else NA_real_,
    max = if (length(x)) max(x) else NA_real_,
    stringsAsFactors = FALSE
  )
}))
write_tsv(gene_count_summary, file.path(analysis_root, "null_gene_count_summary.tsv"))

plot_columns <- c("DESeq2_DEG", "DEXSeq_AS_gene", "DRIMSeq_AS_gene", "SUPPA2_BH_AS_gene")
plot_data <- gene_counts[, plot_columns, drop = FALSE]
if (any(vapply(plot_data, function(z) any(!is.na(z)), logical(1)))) {
  png(
    file.path(analysis_root, "null_gene_count_distributions.png"),
    width = 1800, height = 1200, res = 180
  )
  graphics::par(mar = c(9, 5, 2, 1))
  graphics::boxplot(
    plot_data,
    las = 2,
    ylab = "Significant genes per balanced normoxia split",
    col = c("#4C78A8", "#F58518", "#54A24B", "#E45756"),
    outline = TRUE
  )
  graphics::stripchart(
    plot_data,
    vertical = TRUE,
    method = "jitter",
    add = TRUE,
    pch = 16,
    col = grDevices::adjustcolor("black", alpha.f = 0.45)
  )
  grDevices::dev.off()
}

message("[complete] summaries written to: ", analysis_root)
