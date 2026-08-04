# ============================================================
# Project / Environment
# ============================================================

# Checks if package is loaded
.check_pkg <- function(pkg) {
  for (p in pkg) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Package '", p, "' is required. Load it in the notebook (library(", p, ")).")
    }
  }
  invisible(TRUE)
}

# Returns project root
.get_project_root <- function() {
  root <- Sys.getenv("PROJECT_ROOT")
  if (root == "") {
    if (!requireNamespace("here", quietly = TRUE)) {
      stop("PROJECT_ROOT not set and 'here' not installed.")
    }
    root <- here::here()
  }
  root
}

# Returns results directory
.get_results_dir <- function() {
  results_dir <- Sys.getenv("RESULTS_DIR")

  if (!nzchar(results_dir)) {
    stop(
      "RESULTS_DIR environment variable is not set.\n",
      "Set it in your notebook setup chunk, e.g.:\n",
      "  RESULTS_DIR <- file.path(PROJECT_ROOT, 'results')\n",
      "  Sys.setenv(RESULTS_DIR = RESULTS_DIR)\n"
    )
  }

  if (!dir.exists(results_dir)) {
    stop("RESULTS_DIR does not exist on disk: ", results_dir)
  }

  results_dir
}

# Returns cache directory
.get_cache_dir <- function(resources_dir, subdir = NULL) {
  cache_dir <- file.path(resources_dir, "cache")

  if (!is.null(subdir)) {
    cache_dir <- file.path(cache_dir, subdir)
  }

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_dir
}

# Resolves reference paths from environment
.resolve_ref_from_env <- function(env_file, resources_dir) {
  readRenviron(env_file)

  genome_build <- Sys.getenv("GENOME_BUILD")
  gencode_release <- Sys.getenv("GENCODE_RELEASE")
  genome_gtf <- Sys.getenv("GENOME_GTF")

  if (genome_build == "" || gencode_release == "" || genome_gtf == "") {
    stop("Env file missing one of: GENOME_BUILD, GENCODE_RELEASE, GENOME_GTF")
  }

  genome_dir <- file.path(resources_dir, "genome")
  cache_dir <- file.path(resources_dir, "cache")

  if (!dir.exists(genome_dir)) {
    stop("Genome resources directory not found: ", genome_dir)
  }

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  gtf_path <- file.path(genome_dir, genome_gtf)
  if (!file.exists(gtf_path)) stop("GTF file not found: ", gtf_path)

  gtf_base <- sub("\\.gtf$", "", basename(genome_gtf))
  gtf_base <- sub("\\.gz$", "", gtf_base)

  transcripts_fa <- file.path(genome_dir, paste0(gtf_base, ".transcripts.fa"))
  transcripts_fa_gz <- paste0(transcripts_fa, ".gz")

  fasta_path <- if (file.exists(transcripts_fa)) {
    transcripts_fa
  } else if (file.exists(transcripts_fa_gz)) {
    transcripts_fa_gz
  } else {
    stop(
      "Transcript FASTA not found.\nExpected one of:\n",
      transcripts_fa, "\n",
      transcripts_fa_gz, "\n\n",
      "Run the Salmon preparation script that generates it."
    )
  }

  list(
    resources_dir = resources_dir,
    genome_dir = genome_dir,
    cache_dir = cache_dir,
    genome_build = genome_build,
    gencode_release = gencode_release,
    genome_gtf = genome_gtf,
    gtf_base = gtf_base,
    gtf_path = gtf_path,
    fasta_path = fasta_path
  )
}

# Resolves thread count from environment
.resolve_threads_from_env <- function(env_file) {
  readRenviron(env_file)

  threads <- Sys.getenv("THREADS")

  if (!nzchar(threads)) {
    stop("Env file missing: THREADS")
  }

  threads <- suppressWarnings(as.integer(threads))

  if (is.na(threads) || threads < 1) {
    stop("THREADS must be a positive integer in env file.")
  }

  threads
}

.get_runtime_paths <- function(
  outdir = NULL,
  create_outdir = TRUE
) {
  if (!is.null(outdir) && isTRUE(create_outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }

  list(
    outdir = outdir,
    samples_tsv = file.path(PROJECT_ROOT, "config/metadata/samples.tsv"),
    coldata_tsv = file.path(PROJECT_ROOT, "results/counts/salmon/coldata.tsv"),
    env_file = file.path(PROJECT_ROOT, "config/environment/analysis_environment.env"),
    resources_dir = file.path(PROJECT_ROOT, "resources"),
    indexDir = file.path(PROJECT_ROOT, "resources/salmon_index")
  )
}

.set_outdir <- function(
  outdir = NULL,
  subdir = NULL,
  base_dir = file.path(.get_results_dir(), "plots"),
  create = TRUE
) {
  if (!is.null(outdir) && nzchar(outdir)) {
    resolved <- outdir
  } else {
    if (is.null(subdir) || !nzchar(subdir)) {
      stop("Either `outdir` or `subdir` must be provided.")
    }

    resolved <- file.path(base_dir, subdir)
  }

  resolved <- normalizePath(resolved, mustWork = FALSE)

  if (isTRUE(create)) {
    dir.create(resolved, recursive = TRUE, showWarnings = FALSE)
  }

  resolved
}

# ============================================================
# Utilities
# ===========================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

# Stripts Ensembl version from identifiers
.strip_ens_version <- function(x) {
  # removes anything after first dot: ENSG... .12 -> ENSG...
  sub("\\..*$", "", x)
}

# Read coldata TSV and resolve file paths
.read_coldata_tsv <- function(coldata_tsv,
                              contrast_var = "condition",
                              sample_col = "names",
                              files_col = "files") {
  if (!file.exists(coldata_tsv)) stop("coldata_tsv not found: ", coldata_tsv)

  coldata <- read.delim(coldata_tsv, stringsAsFactors = FALSE, check.names = FALSE)

  required <- c(sample_col, files_col, contrast_var)
  missing_cols <- setdiff(required, colnames(coldata))
  if (length(missing_cols) > 0) {
    stop("coldata_tsv missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (any(!grepl("quant\\.sf$", coldata[[files_col]]))) {
    stop("coldata$", files_col, " must point to Salmon quant.sf files (end with 'quant.sf').")
  }

  project_root <- .get_project_root()

  is_relative <- !grepl("^(/|[A-Za-z]:)", coldata[[files_col]])
  coldata[[files_col]][is_relative] <- file.path(project_root, coldata[[files_col]][is_relative])

  missing_files <- coldata[[files_col]][!file.exists(coldata[[files_col]])]
  if (length(missing_files) > 0) {
    stop(
      "Some quant.sf files do not exist. Examples:\n",
      paste(utils::head(missing_files, 5), collapse = "\n")
    )
  }

  coldata
}

# Pull unique gene IDs from a named list of result data.frames by timepoint
.get_sig_genes_by_tp <- function(sig_by_tp_list) {
  lapply(sig_by_tp_list, function(df) {
    if (is.null(df) || nrow(df) == 0) {
      return(character(0))
    }

    if (!("ensgene" %in% colnames(df))) {
      stop(
        "Expected column 'ensgene' in result table. Available columns: ",
        paste(colnames(df), collapse = ", ")
      )
    }

    unique(stats::na.omit(as.character(df$ensgene)))
  })
}

# Union multiple named timepoint -> gene vectors
.union_gene_lists_by_tp <- function(...) {
  lsts <- list(...)
  tps <- Reduce(union, lapply(lsts, names))

  stats::setNames(lapply(tps, function(tp) {
    unique(unlist(lapply(lsts, function(x) {
      if (tp %in% names(x)) x[[tp]] else character(0)
    }), use.names = FALSE))
  }), tps)
}

# x minus y for named timepoint -> gene vectors
.subtract_gene_lists_by_tp <- function(x, y) {
  tps <- union(names(x), names(y))

  stats::setNames(lapply(tps, function(tp) {
    x_vals <- if (tp %in% names(x)) unique(x[[tp]]) else character(0)
    y_vals <- if (tp %in% names(y)) unique(y[[tp]]) else character(0)
    setdiff(x_vals, y_vals)
  }), tps)
}

# restrict a named timepoint -> gene vector list to a specified set/order of timepoints
.subset_gene_lists_by_tp <- function(x, timepoints) {
  stats::setNames(
    lapply(timepoints, function(tp) {
      if (tp %in% names(x)) x[[tp]] else character(0)
    }),
    timepoints
  )
}

# intersection for named timepoint -> gene vectors
.intersect_gene_lists_by_tp <- function(x, y) {
  tps <- union(names(x), names(y))

  stats::setNames(lapply(tps, function(tp) {
    x_vals <- if (tp %in% names(x)) unique(x[[tp]]) else character(0)
    y_vals <- if (tp %in% names(y)) unique(y[[tp]]) else character(0)
    intersect(x_vals, y_vals)
  }), tps)
}

# Restrict named timepoint -> gene vectors to a named timepoint -> universe vector
.filter_gene_lists_to_universe_by_tp <- function(gene_lists_by_tp, universe_by_tp) {
  tps <- union(names(gene_lists_by_tp), names(universe_by_tp))

  stats::setNames(lapply(tps, function(tp) {
    gene_vals <- if (tp %in% names(gene_lists_by_tp)) unique(gene_lists_by_tp[[tp]]) else character(0)
    universe_vals <- if (tp %in% names(universe_by_tp)) unique(universe_by_tp[[tp]]) else character(0)
    intersect(gene_vals, universe_vals)
  }), tps)
}

# Build timepoint-specific active gene universe from DEseq2
.get_active_gene_universe_by_tp <- function(deseq_results,
                                            timepoints = c("H1", "H3", "H24")) {
  if (is.null(deseq_results$results$tested_universe_by_tp)) {
    stop(
      "deseq_results$results$tested_universe_by_tp not found.\n",
      "Re-run run_deseq() with the updated code that stores tested universes."
    )
  }

  universe_by_tp <- deseq_results$results$tested_universe_by_tp

  missing_tps <- setdiff(timepoints, names(universe_by_tp))
  if (length(missing_tps) > 0) {
    stop(
      "Missing tested universe for timepoint(s): ",
      paste(missing_tps, collapse = ", ")
    )
  }

  universe_by_tp <- universe_by_tp[timepoints]

  universe_by_tp <- lapply(universe_by_tp, function(x) {
    unique(stats::na.omit(.strip_ens_version(as.character(x))))
  })

  universe_by_tp
}

.rbind_fill <- function(x) {
  x <- Filter(Negate(is.null), x)

  if (length(x) == 0) {
    return(data.frame())
  }

  all_cols <- unique(unlist(lapply(x, colnames)))

  x_aligned <- lapply(x, function(df) {
    missing_cols <- setdiff(all_cols, colnames(df))
    if (length(missing_cols) > 0) {
      for (nm in missing_cols) {
        df[[nm]] <- NA
      }
    }
    df[, all_cols, drop = FALSE]
  })

  do.call(rbind, x_aligned)
}

.rbind_fill_df <- function(df_list) {
  df_list <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, df_list)

  if (length(df_list) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }

  all_cols <- unique(unlist(lapply(df_list, colnames)))

  df_list_filled <- lapply(df_list, function(df) {
    missing_cols <- setdiff(all_cols, colnames(df))
    if (length(missing_cols) > 0) {
      for (nm in missing_cols) {
        df[[nm]] <- NA_character_
      }
    }
    df <- df[, all_cols, drop = FALSE]
    rownames(df) <- NULL
    df
  })

  out <- base::do.call(base::rbind, df_list_filled)
  rownames(out) <- NULL
  out
}

# ============================================================
# Annotations
# ============================================================

# Checks for linkedTXome (tximeta)
.ensure_linked_txome <- function(indexDir,
                                 env_file,
                                 resources_dir,
                                 source = "LocalGENCODE",
                                 organism = "Homo sapiens") {
  .check_pkg("tximeta")

  if (!dir.exists(indexDir)) stop("Salmon index directory not found: ", indexDir)

  ref <- .resolve_ref_from_env(env_file, resources_dir)

  linked_json <- file.path(dirname(indexDir), paste0(basename(indexDir), ".json"))
  if (file.exists(linked_json)) {
    message("linkedTxome already exists; skipping makeLinkedTxome(): ", linked_json)
    return(invisible(c(ref, list(indexDir = indexDir, linked_json = linked_json, created = FALSE))))
  }

  message("Creating linkedTxome...: ", linked_json)
  message("Genome: ", ref$genome_build)
  message("GENCODE Release: ", ref$gencode_release)
  message("Transcript FASTA: ", ref$fasta_path)
  message("GTF: ", ref$gtf_path)

  tximeta::makeLinkedTxome(
    indexDir  = indexDir,
    source    = source,
    organism  = organism,
    release   = ref$gencode_release,
    genome    = ref$genome_build,
    fasta     = ref$fasta_path,
    gtf       = ref$gtf_path
  )

  invisible(c(ref, list(indexDir = indexDir, linked_json = linked_json, created = TRUE)))
}

# Get a simple annotation table from annotables
.get_annot <- function() {
  .check_pkg("annotables")
  annot_df <- annotables::grch38[, c("ensgene", "symbol", "description")]
  annot_df <- annot_df[!duplicated(annot_df$ensgene), , drop = FALSE]
  annot_df
}

# Attaches annoation to results
.attach_gene_annot <- function(out, annot_df = NULL) {
  if (is.null(annot_df)) {
    annot_df <- .get_annot()
  }

  idx <- match(out$ensgene, annot_df$ensgene)
  cbind(
    out,
    annot_df[idx, setdiff(colnames(annot_df), "ensgene"), drop = FALSE]
  )
}

# Annotate Ensembl gene IDs using annotables
.annotate_ensgenes <- function(ensgenes, annot_df = NULL) {
  if (is.null(annot_df)) {
    annot_df <- .get_annot()
  }

  out <- data.frame(
    ensgene = unique(.strip_ens_version(ensgenes)),
    stringsAsFactors = FALSE
  )

  idx <- match(out$ensgene, annot_df$ensgene)
  out <- cbind(
    out,
    annot_df[idx, setdiff(colnames(annot_df), "ensgene"), drop = FALSE]
  )

  out
}

# Get or create transcript-level tximeta SE
.get_se <- function(coldata_tsv,
                    env_file,
                    resources_dir,
                    indexDir,
                    force = FALSE,
                    skipMeta = FALSE) {
  .check_pkg("tximeta")
  .check_pkg("SummarizedExperiment")
  .check_pkg("S4Vectors")

  cache_dir <- .get_cache_dir(resources_dir)

  se_rds <- file.path(cache_dir, "tximeta_se.rds")
  if (file.exists(se_rds) && !force) {
    message("Loading cached transcript-level SummarizedExperiment: ", se_rds)
    return(readRDS(se_rds))
  }

  .ensure_linked_txome(
    indexDir = indexDir,
    env_file = env_file,
    resources_dir = resources_dir
  )

  coldata <- .read_coldata_tsv(
    coldata_tsv = coldata_tsv,
    contrast_var = "condition",
    sample_col = "names",
    files_col = "files"
  )

  se <- suppressWarnings(
    tximeta::tximeta(coldata, skipMeta = skipMeta)
  )

  # Convert character colData columns to factors so downstream design formulas behave cleanly
  cd <- as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)
  is_char <- vapply(cd, is.character, logical(1))
  if (any(is_char)) {
    for (nm in names(cd)[is_char]) {
      cd[[nm]] <- factor(cd[[nm]])
    }
    SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(cd)
  }

  saveRDS(se, se_rds)
  message("Saved transcript-level SummarizedExperiment to: ", se_rds)
  se
}

# Get or create gene-summarized tximeta object
.get_gse <- function(coldata_tsv,
                     env_file,
                     resources_dir,
                     indexDir,
                     force = FALSE,
                     skipMeta = FALSE,
                     assignRanges = "abundant") {
  .check_pkg("tximeta")
  .check_pkg("SummarizedExperiment")
  .check_pkg("S4Vectors")

  cache_dir <- .get_cache_dir(resources_dir)
  gse_rds <- file.path(cache_dir, "tximeta_gse.rds")

  if (file.exists(gse_rds) && !force) {
    message("Loading cached gene-level SummarizedExperiment: ", gse_rds)
    return(readRDS(gse_rds))
  }

  se <- .get_se(
    coldata_tsv = coldata_tsv,
    env_file = env_file,
    resources_dir = resources_dir,
    indexDir = indexDir,
    force = force,
    skipMeta = skipMeta
  )

  gse <- tximeta::summarizeToGene(se, assignRanges = assignRanges)

  # Convert character colData columns to factors so downstream design formulas behave cleanly
  cd <- as.data.frame(SummarizedExperiment::colData(gse), stringsAsFactors = FALSE)
  is_char <- vapply(cd, is.character, logical(1))
  if (any(is_char)) {
    for (nm in names(cd)[is_char]) {
      cd[[nm]] <- factor(cd[[nm]])
    }
    SummarizedExperiment::colData(gse) <- S4Vectors::DataFrame(cd)
  }

  saveRDS(gse, gse_rds)
  message("Saved gene-level SummarizedExperiment to: ", gse_rds)
  gse
}

.get_txdb <- function(env_file,
                      resources_dir,
                      overwrite = FALSE) {
  .check_pkg("txdbmaker")
  .check_pkg("GenomeInfoDb")
  .check_pkg("AnnotationDbi")

  ref <- .resolve_ref_from_env(env_file, resources_dir)
  genome_dir <- ref$genome_dir
  txdb_sqlite <- file.path(genome_dir, paste0(ref$gtf_base, ".txdb.sqlite"))

  if (file.exists(txdb_sqlite) && !overwrite) {
    message("Loading cached TxDb: ", txdb_sqlite)
    txdb <- AnnotationDbi::loadDb(txdb_sqlite)
    return(invisible(c(
      ref,
      list(
        txdb = txdb,
        txdb_sqlite = txdb_sqlite,
        created = FALSE
      )
    )))
  }

  message("Building TxDb from GTF (this happens once)...")
  txdb <- suppressWarnings(txdbmaker::makeTxDbFromGFF(ref$gtf_path))
  GenomeInfoDb::genome(GenomeInfoDb::seqinfo(txdb)) <- ref$genome_build

  AnnotationDbi::saveDb(txdb, file = txdb_sqlite)
  message("TxDb saved to: ", txdb_sqlite)

  invisible(c(
    ref,
    list(
      txdb = txdb,
      txdb_sqlite = txdb_sqlite,
      created = TRUE
    )
  ))
}

# Get or create tx2gene mapping from the cached/imported SE annotation
.get_tx2gene <- function(env_file,
                         resources_dir,
                         force = FALSE) {
  .check_pkg("AnnotationDbi")

  cache_dir <- .get_cache_dir(resources_dir)
  tx2gene_rds <- file.path(cache_dir, "tx2gene.rds")
  tx2gene_tsv <- file.path(cache_dir, "tx2gene.tsv")

  if (file.exists(tx2gene_rds) && !force) {
    message("Loading cached tx2gene: ", tx2gene_rds)
    return(readRDS(tx2gene_rds))
  }

  txdb_info <- .get_txdb(
    env_file = env_file,
    resources_dir = resources_dir,
    overwrite = force
  )

  txdb <- txdb_info$txdb

  k <- AnnotationDbi::keys(txdb, keytype = "TXNAME")
  tx2gene <- AnnotationDbi::select(
    txdb,
    keys = k,
    columns = "GENEID",
    keytype = "TXNAME"
  )

  tx2gene <- tx2gene[!is.na(tx2gene$TXNAME) & !is.na(tx2gene$GENEID), c("TXNAME", "GENEID")]
  tx2gene <- tx2gene[!duplicated(tx2gene$TXNAME), , drop = FALSE]
  rownames(tx2gene) <- NULL

  saveRDS(tx2gene, tx2gene_rds)
  utils::write.table(
    tx2gene,
    file = tx2gene_tsv,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  message("Saved tx2gene to: ", tx2gene_rds)
  tx2gene
}

# ============================================================
# GO Enrichment
# ============================================================

# Retrieves OBO file matching PANTHER release
.get_go_obo <- function(resources_dir,
                        dataset_id = "GO:0008150",
                        fallback_date = "2025-10-10",
                        force = FALSE) {
  .check_pkg("jsonlite")
  .check_pkg("xml2")

  go_dir <- file.path(resources_dir, "ontology")
  dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)

  # get PANTHER GO release metadata
  meta <- tryCatch(
    {
      x <- jsonlite::fromJSON(
        "https://pantherdb.org/services/oai/pantherdb/supportedannotdatasets"
      )

      adt <- as.data.frame(
        x$search$annotation_data_sets$annotation_data_type,
        stringsAsFactors = FALSE
      )

      adt[adt$id == dataset_id, c("id", "label", "release_date", "version"), drop = FALSE]
    },
    error = function(e) NULL
  )

  if (!is.null(meta) && nrow(meta) > 0 && nzchar(meta$release_date[1])) {
    target_date <- as.Date(meta$release_date[1])
    dataset_label <- meta$label[1]
    version <- meta$version[1]
  } else {
    target_date <- as.Date(fallback_date)
    dataset_label <- NA_character_
    version <- NA_character_
  }

  # find closest GO archive release
  matched_release <- tryCatch(
    {
      archive_url <- "https://release.geneontology.org/"
      page <- xml2::read_html(archive_url)

      hrefs <- xml2::xml_attr(xml2::xml_find_all(page, ".//a"), "href")
      hrefs <- hrefs[!is.na(hrefs)]

      rels <- unique(gsub("/$", "", hrefs[grepl("^\\d{4}-\\d{2}-\\d{2}/?$", hrefs)]))
      rel_dates <- as.Date(rels)

      rel_dates <- sort(unique(rel_dates[!is.na(rel_dates)]))

      if (length(rel_dates) == 0 || is.na(target_date)) {
        as.Date(fallback_date)
      } else {
        after_idx <- which(rel_dates >= target_date)

        if (length(after_idx) > 0) {
          rel_dates[min(after_idx)]
        } else {
          before_idx <- which(rel_dates < target_date)
          if (length(before_idx) > 0) rel_dates[max(before_idx)] else as.Date(fallback_date)
        }
      }
    },
    error = function(e) {
      as.Date(fallback_date)
    }
  )

  # resolve URL + cache location
  matched_release <- as.character(matched_release)

  obo_url <- paste0(
    "https://release.geneontology.org/",
    matched_release,
    "/ontology/go-basic.obo"
  )

  obo_file <- file.path(
    go_dir,
    paste0("go-basic_", matched_release, ".obo")
  )

  # download file (if needed)

  if (!file.exists(obo_file) || isTRUE(force)) {
    message("[PANTHER] Downloading: ", obo_url)

    utils::download.file(
      url = obo_url,
      destfile = obo_file,
      mode = "wb",
      quiet = FALSE
    )
  } else {
    message("[PANTHER] Using cached file: ", basename(obo_file))
  }

  # return metadata
  list(
    dataset_id = dataset_id,
    dataset_label = dataset_label,
    panther_release_date = as.character(target_date),
    matched_go_release_date = matched_release,
    version = version,
    obo_url = obo_url,
    obo_file = obo_file
  )
}

# Parses PANTHER mapping file to generate relationship between GO terms and gene symbols
.parse_panther_gene_annotations <- function(gene_obj) {
  if (is.null(gene_obj) || !is.list(gene_obj)) {
    return(data.frame(stringsAsFactors = FALSE))
  }

  if (!("annotation_type_list" %in% names(gene_obj))) {
    return(data.frame(stringsAsFactors = FALSE))
  }

  symbol <- if ("mapped_id_list" %in% names(gene_obj)) {
    as.character(gene_obj$mapped_id_list)
  } else {
    NA_character_
  }

  adt <- gene_obj$annotation_type_list$annotation_data_type
  if (is.null(adt) || !is.list(adt)) {
    return(data.frame(stringsAsFactors = FALSE))
  }

  # single block -> list of blocks
  if (!is.null(names(adt)) && any(names(adt) %in% c("content", "annotation_list", "release_version"))) {
    adt <- list(adt)
  }

  out <- lapply(adt, function(block) {
    if (is.null(block) || !is.list(block)) {
      return(NULL)
    }

    content <- if ("content" %in% names(block)) as.character(block$content) else NA_character_
    release_version <- if ("release_version" %in% names(block)) as.character(block$release_version) else NA_character_

    anns <- block$annotation_list$annotation
    if (is.null(anns) || !is.list(anns)) {
      return(NULL)
    }

    # single annotation -> list of annotations
    if (!is.null(names(anns)) && any(names(anns) %in% c("id", "name"))) {
      anns <- list(anns)
    }

    ann_rows <- lapply(anns, function(ann) {
      if (is.null(ann) || !is.list(ann)) {
        return(NULL)
      }

      data.frame(
        symbol = symbol,
        annot_dataset = content,
        release_version = release_version,
        pathway_id = if ("id" %in% names(ann)) as.character(ann$id) else NA_character_,
        pathway_label = if ("name" %in% names(ann)) as.character(ann$name) else NA_character_,
        stringsAsFactors = FALSE
      )
    })

    .rbind_fill_df(ann_rows)
  })

  .rbind_fill_df(out)
}

# Retrieves PANTHER mapping object for given GO datasets
.get_panther_mapping <- function(res,
                                 annot_datasets = c("GO:0008150", "GO:0003674"),
                                 annot_df = NULL,
                                 organism = 9606,
                                 resources_dir,
                                 force = FALSE,
                                 save_tables = FALSE) {
  .check_pkg("rbioapi")

  cache_dir <- .get_cache_dir(resources_dir)
  cache_file <- file.path(cache_dir, "term2gene_mapping.rds")

  if (file.exists(cache_file) && !force) {
    message("[PANTHER] Loading cached mapping: ", basename(cache_file))
    return(readRDS(cache_file))
  }

  if (is.null(annot_df)) {
    annot_df <- res$annotation$annot_df
    if (is.null(annot_df)) {
      annot_df <- .get_annot()
    }
  }

  tested_universe_by_tp <- res$gene_sets$tested_universe_by_tp
  if (is.null(tested_universe_by_tp) || length(tested_universe_by_tp) == 0) {
    stop("res$gene_sets$tested_universe_by_tp is empty.")
  }

  universe_all_ens <- unique(unlist(tested_universe_by_tp, use.names = FALSE))
  universe_all_ens <- unique(stats::na.omit(as.character(universe_all_ens)))

  universe_all_symbols <- unique(stats::na.omit(as.character(
    annot_df$symbol[match(universe_all_ens, annot_df$ensgene)]
  )))
  universe_all_symbols <- universe_all_symbols[universe_all_symbols != ""]

  chunks <- if (length(universe_all_symbols) == 0) {
    list(character(0))
  } else {
    split(universe_all_symbols, ceiling(seq_along(universe_all_symbols) / 1000))
  }

  message("[PANTHER] Universe size (Ensembl): ", length(universe_all_ens))
  message("[PANTHER] Universe size (symbols): ", length(universe_all_symbols))
  message("[PANTHER] Mapping ", length(universe_all_symbols), " symbols in ", length(chunks), " chunk(s)")

  parsed_chunks <- lapply(seq_along(chunks), function(i) {
    chunk_syms <- chunks[[i]]
    message("[PANTHER mapping]   chunk ", i, "/", length(chunks), " (", length(chunk_syms), " genes)")

    tmp <- rbioapi::rba_panther_mapping(
      genes = chunk_syms,
      organism = organism
    )

    if (is.null(tmp$mapped_genes) || is.null(tmp$mapped_genes$gene) || length(tmp$mapped_genes$gene) == 0) {
      return(data.frame(stringsAsFactors = FALSE))
    }

    gene_entries <- tmp$mapped_genes$gene

    if (!is.list(gene_entries)) {
      return(data.frame(stringsAsFactors = FALSE))
    }

    if (!is.null(names(gene_entries)) &&
      any(names(gene_entries) %in% c("mapped_id_list", "accession", "annotation_type_list"))) {
      gene_entries <- list(gene_entries)
    }

    parsed_list <- lapply(gene_entries, .parse_panther_gene_annotations)
    parsed_list <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, parsed_list)

    .rbind_fill_df(parsed_list)
  })

  mapping_long <- .rbind_fill_df(parsed_chunks)

  if (nrow(mapping_long) > 0) {
    mapping_long <- mapping_long[
      mapping_long$annot_dataset %in% annot_datasets &
        !is.na(mapping_long$symbol) & mapping_long$symbol != "" &
        !is.na(mapping_long$pathway_id) & mapping_long$pathway_id != "", ,
      drop = FALSE
    ]
    mapping_long <- unique(mapping_long)
  }

  mapping_summary <- if (nrow(mapping_long) == 0) {
    data.frame()
  } else {
    tmp <- aggregate(
      symbol ~ annot_dataset + pathway_id + pathway_label,
      data = mapping_long,
      FUN = function(x) list(sort(unique(as.character(x))))
    )
    colnames(tmp)[colnames(tmp) == "symbol"] <- "symbols"
    tmp$pathway_size <- vapply(tmp$symbols, length, integer(1))
    tmp
  }

  out <- list(
    meta = list(
      annot_datasets = annot_datasets,
      organism = organism
    ),
    universe = list(
      ensgene = universe_all_ens,
      symbol = universe_all_symbols
    ),
    long = mapping_long,
    summary = mapping_summary
  )

  if (isTRUE(save_tables)) {
    if (nrow(mapping_long) > 0) {
      utils::write.table(
        mapping_long,
        file = file.path(cache_dir, "panther_mapping_long.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
    }

    if (nrow(mapping_summary) > 0) {
      tmp <- mapping_summary
      tmp$symbols <- vapply(tmp$symbols, paste, collapse = ",", FUN.VALUE = character(1))
      utils::write.table(
        tmp,
        file = file.path(cache_dir, "panther_mapping_summary.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
    }
  }

  saveRDS(out, cache_file)
  out
}

# Combines OBO and PANTHER mapping to propagate through the ontology hierarchy and generate a comprehensive term2gene relationship
.get_term2gene <- function(mapping_res,
                           resources_dir,
                           fallback_go_release_date = "2025-10-10",
                           relationships = c("is_a", "part_of"),
                           force_obo = FALSE,
                           force = FALSE,
                           save_tables = FALSE) {
  .check_pkg("ontologyIndex")

  cache_dir <- .get_cache_dir(resources_dir)
  cache_file <- file.path(cache_dir, "term2gene_propagated.rds")

  if (file.exists(cache_file) && !force) {
    message("[PANTHER] Loading cached propagated term2gene: ", basename(cache_file))
    return(readRDS(cache_file))
  }

  mapping_long <- mapping_res$long
  annot_datasets <- unique(as.character(mapping_long$annot_dataset))

  term2gene_by_dataset <- stats::setNames(vector("list", length(annot_datasets)), annot_datasets)
  obo_info_by_dataset <- stats::setNames(vector("list", length(annot_datasets)), annot_datasets)

  for (ds in annot_datasets) {
    message("[PANTHER] Dataset: ", ds)

    obo_info <- .get_go_obo(
      resources_dir = resources_dir,
      dataset_id = ds,
      fallback_date = fallback_go_release_date,
      force = force_obo
    )
    obo_info_by_dataset[[ds]] <- obo_info

    go <- ontologyIndex::get_ontology(
      obo_info$obo_file,
      propagate_relationships = relationships
    )

    df <- mapping_long[mapping_long$annot_dataset == ds, c("pathway_id", "pathway_label", "symbol"), drop = FALSE]
    df <- unique(df)
    df <- df[!is.na(df$pathway_id) & df$pathway_id != "", , drop = FALSE]
    df <- df[!is.na(df$symbol) & df$symbol != "", , drop = FALSE]
    df <- df[df$pathway_id %in% go$id, , drop = FALSE]

    if (nrow(df) == 0) {
      term2gene_by_dataset[[ds]] <- list(long = data.frame(), summary = data.frame())
      next
    }

    direct_label_map <- stats::setNames(df$pathway_label, df$pathway_id)
    direct_label_map <- direct_label_map[!duplicated(names(direct_label_map))]

    propagated_list <- lapply(seq_len(nrow(df)), function(i) {
      term <- df$pathway_id[i]
      sym <- df$symbol[i]

      anc <- go$ancestors[[term]]
      if (is.null(anc)) anc <- character(0)

      expanded_terms <- unique(c(term, anc))
      expanded_terms <- expanded_terms[expanded_terms %in% go$id]

      data.frame(
        pathway_id = expanded_terms,
        symbol = sym,
        stringsAsFactors = FALSE
      )
    })

    term2gene_long <- .rbind_fill_df(propagated_list)
    term2gene_long <- unique(term2gene_long)

    go_name_map <- go$name
    term2gene_long$pathway_label <- unname(go_name_map[term2gene_long$pathway_id])

    missing_lab <- is.na(term2gene_long$pathway_label) | term2gene_long$pathway_label == ""
    if (any(missing_lab)) {
      term2gene_long$pathway_label[missing_lab] <- unname(
        direct_label_map[term2gene_long$pathway_id[missing_lab]]
      )
    }

    term2gene_long <- term2gene_long[, c("pathway_id", "pathway_label", "symbol"), drop = FALSE]
    term2gene_long <- unique(term2gene_long)
    term2gene_long$annot_dataset <- ds

    term2gene_summary <- aggregate(
      symbol ~ annot_dataset + pathway_id + pathway_label,
      data = term2gene_long,
      FUN = function(x) list(sort(unique(as.character(x))))
    )
    colnames(term2gene_summary)[colnames(term2gene_summary) == "symbol"] <- "symbols"
    term2gene_summary$pathway_size <- vapply(term2gene_summary$symbols, length, integer(1))

    term2gene_by_dataset[[ds]] <- list(
      long = term2gene_long,
      summary = term2gene_summary
    )
  }

  term2gene_long_all <- .rbind_fill_df(lapply(term2gene_by_dataset, `[[`, "long"))
  term2gene_summary_all <- .rbind_fill_df(lapply(term2gene_by_dataset, `[[`, "summary"))

  out <- list(
    meta = list(
      relationships = relationships
    ),
    ontology = obo_info_by_dataset,
    long = term2gene_long_all,
    summary = term2gene_summary_all
  )

  if (isTRUE(save_tables)) {
    if (nrow(term2gene_long_all) > 0) {
      utils::write.table(
        term2gene_long_all,
        file = file.path(cache_dir, "term2gene_long.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
    }

    if (nrow(term2gene_summary_all) > 0) {
      tmp <- term2gene_summary_all
      tmp$symbols <- vapply(tmp$symbols, paste, collapse = ",", FUN.VALUE = character(1))
      utils::write.table(
        tmp,
        file = file.path(cache_dir, "term2gene_summary.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
    }
  }

  saveRDS(out, cache_file)
  out
}

# Run PANTHER GO BP enrichment from Ensembl IDs using annotables symbols
.run_panther_enrich_from_ensgenes <- function(query_ensgenes,
                                              universe_ensgenes,
                                              annot_df = NULL,
                                              organism = 9606,
                                              ref_organism = 9606,
                                              annot_dataset = "GO:0008150",
                                              cutoff = 0.05,
                                              cache_file = NULL,
                                              force = FALSE) {
  .check_pkg("rbioapi")

  if (!is.null(cache_file) && file.exists(cache_file) && !force) {
    message("[PANTHER enrich] Loading cached result: ", basename(cache_file))
    return(readRDS(cache_file))
  }

  if (is.null(annot_df)) {
    annot_df <- .get_annot()
  }

  query_annot <- .annotate_ensgenes(query_ensgenes, annot_df = annot_df)
  universe_annot <- .annotate_ensgenes(universe_ensgenes, annot_df = annot_df)

  query_symbols <- unique(stats::na.omit(as.character(query_annot$symbol)))
  universe_symbols <- unique(stats::na.omit(as.character(universe_annot$symbol)))

  if (length(query_symbols) == 0) {
    out <- list(
      query_annot = query_annot,
      universe_annot = universe_annot,
      result = data.frame(),
      raw = NULL
    )

    if (!is.null(cache_file)) {
      dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
      saveRDS(out, cache_file)
    }

    return(out)
  }

  enr <- rbioapi::rba_panther_enrich(
    genes = query_symbols,
    ref_genes = universe_symbols,
    organism = organism,
    ref_organism = ref_organism,
    annot_dataset = annot_dataset,
    cutoff = cutoff
  )

  res <- enr[["result"]]
  if (is.null(res)) {
    res <- data.frame()
  } else if (nrow(res) > 0) {
    res$gene_ratio <- res$number_in_list / res$number_in_reference
  }

  out <- list(
    query_annot = query_annot,
    universe_annot = universe_annot,
    result = res,
    raw = enr
  )

  if (!is.null(cache_file)) {
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, cache_file)
  }

  out
}

# ============================================================
# Plotting
# ============================================================

plot_feature_overlap_across_timepoints <- function(
  features_by_tp,
  title = NULL,
  shade_label = "TX"
) {
  .check_pkg(c(
    "eulerr",
    "grid"
  ))

  timepoints <- c("H1", "H3", "H24")

  features_by_tp <- lapply(timepoints, function(tp) {
    unique(stats::na.omit(features_by_tp[[tp]]))
  })
  names(features_by_tp) <- timepoints

  h1 <- features_by_tp[["H1"]]
  h3 <- features_by_tp[["H3"]]
  h24 <- features_by_tp[["H24"]]

  tp_shaded_colors <- .get_timepoint_shaded_colors(
    timepoints = timepoints,
    labels = c("TX", "OVERLAP", "GENE")
  )

  tp_colors <- vapply(
    tp_shaded_colors,
    function(x) unname(x[[shade_label]]),
    character(1)
  )

  fit <- eulerr::euler(
    c(
      "1H" = length(h1),
      "3H" = length(h3),
      "24H" = length(h24),
      "1H&3H" = length(intersect(h1, h3)),
      "1H&24H" = length(intersect(h1, h24)),
      "3H&24H" = length(intersect(h3, h24)),
      "1H&3H&24H" = length(Reduce(intersect, list(h1, h3, h24)))
    ),
    input = "union"
  )

  p <- plot(
    fit,
    fills = unname(tp_colors[timepoints]),
    quantities = TRUE,
    legend = TRUE,
    main = title
  )

  invisible(list(
    plot = p,
    fit = fit,
    features_by_tp = features_by_tp,
    counts = vapply(features_by_tp, length, integer(1))
  ))
}
