#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  required_pkgs <- c("jsonlite", "tximport", "DESeq2", "rhdf5")
  missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(paste("missing required R packages:", paste(missing, collapse = ", ")), call. = FALSE)
  }
}))

parse_args <- function(raw_args) {
  parsed <- list(contrast = character())
  i <- 1L
  while (i <= length(raw_args)) {
    key <- raw_args[[i]]
    if (!startsWith(key, "--")) {
      stop(paste("unexpected positional argument:", key), call. = FALSE)
    }
    if (i == length(raw_args)) {
      stop(paste("missing value for", key), call. = FALSE)
    }
    value <- raw_args[[i + 1L]]
    name <- substring(key, 3L)
    if (name == "contrast") {
      parsed$contrast <- c(parsed$contrast, value)
    } else {
      parsed[[name]] <- value
    }
    i <- i + 2L
  }
  parsed
}

required_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || identical(value, "")) {
    stop(paste("missing required argument --", name, sep = ""), call. = FALSE)
  }
  value
}

parse_bool <- function(value, default = FALSE) {
  if (is.null(value) || identical(value, "")) {
    return(default)
  }
  normalized <- tolower(trimws(value))
  if (normalized %in% c("1", "true", "t", "yes", "y")) {
    return(TRUE)
  }
  if (normalized %in% c("0", "false", "f", "no", "n")) {
    return(FALSE)
  }
  stop(paste("invalid boolean value:", value), call. = FALSE)
}

split_contrast <- function(raw) {
  parts <- strsplit(raw, "|", fixed = TRUE)[[1]]
  if (length(parts) < 3L) {
    stop(paste("invalid contrast spec:", raw), call. = FALSE)
  }
  if (length(parts) < 4L || identical(parts[[4]], "")) {
    parts[[4]] <- paste(parts[[1]], parts[[2]], "vs", parts[[3]], sep = "_")
  }
  list(factor = parts[[1]], level_a = parts[[2]], level_b = parts[[3]], name = parts[[4]])
}

design_variables <- function(formula_string) {
  terms <- attr(stats::terms(stats::as.formula(formula_string)), "term.labels")
  if (length(terms) == 0L) {
    return(character())
  }
  unique(trimws(unlist(strsplit(terms, "[:*+]"))))
}

sanitize_name <- function(value) {
  gsub("[^A-Za-z0-9_.-]", "_", value)
}

normalize_stable_id <- function(value) {
  sub("\\.[0-9]+$", "", value)
}

load_gene_annotation <- function(path) {
  if (is.null(path) || identical(path, "")) {
    return(NULL)
  }
  annotation <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"gene_id" %in% colnames(annotation)) {
    stop("gene annotation must contain gene_id column", call. = FALSE)
  }
  if (!"gene_name" %in% colnames(annotation)) {
    annotation$gene_name <- ""
  }
  if (!"gene_biotype" %in% colnames(annotation)) {
    annotation$gene_biotype <- ""
  }
  annotation$gene_id <- normalize_stable_id(annotation$gene_id)
  annotation <- annotation[, c("gene_id", "gene_name", "gene_biotype")]
  annotation <- annotation[!duplicated(annotation$gene_id), , drop = FALSE]
  rownames(annotation) <- annotation$gene_id
  annotation
}

annotate_table <- function(df, annotation) {
  if (is.null(annotation) || !"gene_id" %in% colnames(df)) {
    return(df)
  }
  idx <- match(df$gene_id, annotation$gene_id)
  gene_name <- annotation$gene_name[idx]
  gene_biotype <- annotation$gene_biotype[idx]
  gene_name[is.na(gene_name)] <- ""
  gene_biotype[is.na(gene_biotype)] <- ""
  payload_cols <- setdiff(colnames(df), "gene_id")
  cbind(
    gene_id = df$gene_id,
    gene_symbol = gene_name,
    gene_biotype = gene_biotype,
    df[, payload_cols, drop = FALSE]
  )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  samplesheet_path <- required_arg(args, "samplesheet")
  quant_root <- required_arg(args, "quant-root")
  engine <- required_arg(args, "engine")
  tx2gene_path <- required_arg(args, "tx2gene")
  gene_annotation_path <- args[["gene-annotation"]]
  design_raw <- required_arg(args, "design")
  transform <- required_arg(args, "transform")
  counts_from_abundance <- required_arg(args, "counts-from-abundance")
  ignore_tx_version <- parse_bool(args[["ignore-tx-version"]], default = TRUE)
  ignore_after_bar <- parse_bool(args[["ignore-after-bar"]], default = FALSE)
  outdir <- required_arg(args, "outdir")
  manifest_out <- required_arg(args, "manifest-out")

  started_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)

  samplesheet <- utils::read.delim(samplesheet_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"run_accession" %in% colnames(samplesheet)) {
    stop("samplesheet must contain run_accession column", call. = FALSE)
  }

  keep_idx <- integer()
  files <- character()
  for (i in seq_len(nrow(samplesheet))) {
    run_id <- samplesheet$run_accession[[i]]
    abundance <- file.path(quant_root, run_id, engine, "abundance.h5")
    if (file.exists(abundance) && file.info(abundance)$size > 0) {
      keep_idx <- c(keep_idx, i)
      files <- c(files, abundance)
    }
  }
  if (length(files) < 2L) {
    stop("need at least 2 quantified runs for DESeq2 stage", call. = FALSE)
  }
  samplesheet <- samplesheet[keep_idx, , drop = FALSE]
  rownames(samplesheet) <- samplesheet$run_accession
  names(files) <- samplesheet$run_accession
  for (column in design_variables(design_raw)) {
    if (!column %in% colnames(samplesheet)) {
      next
    }
    if (is.character(samplesheet[[column]])) {
      samplesheet[[column]] <- as.factor(samplesheet[[column]])
    }
  }

  tx2gene <- utils::read.delim(tx2gene_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(tx2gene) < 2L) {
    stop("tx2gene must contain at least two columns", call. = FALSE)
  }
  tx2gene <- tx2gene[, 1:2]
  colnames(tx2gene) <- c("TXNAME", "GENEID")
  tx2gene$TXNAME <- normalize_stable_id(tx2gene$TXNAME)
  tx2gene$GENEID <- normalize_stable_id(tx2gene$GENEID)
  tx2gene <- tx2gene[!duplicated(tx2gene$TXNAME), , drop = FALSE]
  gene_annotation <- load_gene_annotation(gene_annotation_path)

  tximport_stdout <- utils::capture.output({
    txi <- suppressMessages(suppressWarnings(tximport::tximport(
      files,
      type = "kallisto",
      tx2gene = tx2gene,
      countsFromAbundance = counts_from_abundance,
      ignoreTxVersion = ignore_tx_version,
      ignoreAfterBar = ignore_after_bar
    )))
  })

  design_formula <- stats::as.formula(design_raw)
  dds <- suppressMessages(suppressWarnings(DESeq2::DESeqDataSetFromTximport(
    txi = txi,
    colData = samplesheet,
    design = design_formula
  )))
  dds <- suppressMessages(DESeq2::DESeq(dds, quiet = TRUE))

  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }

  gene_counts <- DESeq2::counts(dds, normalized = FALSE)
  gene_norm <- DESeq2::counts(dds, normalized = TRUE)
  gene_counts_df <- cbind(gene_id = rownames(gene_counts), as.data.frame(gene_counts))
  utils::write.table(
    gene_counts_df,
    file = file.path(outdir, "gene_counts.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    annotate_table(gene_counts_df, gene_annotation),
    file = file.path(outdir, "gene_counts_annotated.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  gene_norm_df <- cbind(gene_id = rownames(gene_norm), as.data.frame(gene_norm))
  utils::write.table(
    gene_norm_df,
    file = file.path(outdir, "gene_norm_counts.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    annotate_table(gene_norm_df, gene_annotation),
    file = file.path(outdir, "gene_norm_counts_annotated.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  if (identical(transform, "rlog")) {
    transformed <- DESeq2::rlog(dds, blind = TRUE)
  } else {
    transformed <- DESeq2::vst(dds, blind = TRUE)
  }
  vst_mat <- SummarizedExperiment::assay(transformed)
  vst_df <- cbind(gene_id = rownames(vst_mat), as.data.frame(vst_mat))
  utils::write.table(
    vst_df,
    file = file.path(outdir, "vst.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    annotate_table(vst_df, gene_annotation),
    file = file.path(outdir, "vst_annotated.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  saveRDS(transformed, file = file.path(outdir, "vst.rds"))

  size_factors <- DESeq2::sizeFactors(dds)
  utils::write.table(
    data.frame(run_accession = names(size_factors), size_factor = as.numeric(size_factors)),
    file = file.path(outdir, "size_factors.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  contrast_specs <- lapply(args$contrast, split_contrast)
  de_outputs <- character()
  for (contrast in contrast_specs) {
    res <- DESeq2::results(
      dds,
      contrast = c(contrast$factor, contrast$level_a, contrast$level_b)
    )
    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)
    out_name <- paste0("de_", sanitize_name(contrast$name), ".tsv")
    out_path <- file.path(outdir, out_name)
    utils::write.table(
      res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))],
      file = out_path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    utils::write.table(
      annotate_table(res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))], gene_annotation),
      file = file.path(outdir, paste0("de_", sanitize_name(contrast$name), "_annotated.tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    de_outputs <- c(de_outputs, out_path)
  }

  manifest <- list(
    schema_version = 1L,
    stage = "deseq2",
    id = basename(outdir),
    started_at = started_at,
    finished_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    parameters = list(
      design = design_raw,
      transform = transform,
      counts_from_abundance = counts_from_abundance,
      ignore_tx_version = ignore_tx_version,
      ignore_after_bar = ignore_after_bar,
      contrasts = args$contrast,
      samples_used = nrow(samplesheet)
    ),
    tool_versions = list(
      R = R.version.string,
      DESeq2 = as.character(utils::packageVersion("DESeq2")),
      tximport = as.character(utils::packageVersion("tximport"))
    ),
    outputs = list(
      gene_counts = file.path(outdir, "gene_counts.tsv"),
      gene_counts_annotated = file.path(outdir, "gene_counts_annotated.tsv"),
      gene_norm_counts = file.path(outdir, "gene_norm_counts.tsv"),
      gene_norm_counts_annotated = file.path(outdir, "gene_norm_counts_annotated.tsv"),
      vst = file.path(outdir, "vst.tsv"),
      vst_annotated = file.path(outdir, "vst_annotated.tsv"),
      vst_rds = file.path(outdir, "vst.rds"),
      size_factors = file.path(outdir, "size_factors.tsv"),
      de_tables = de_outputs
    ),
    captured_stdout = tximport_stdout,
    session_info = utils::capture.output(sessionInfo())
  )
  dir.create(dirname(manifest_out), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(manifest, path = manifest_out, pretty = TRUE, auto_unbox = TRUE)
}

main()
