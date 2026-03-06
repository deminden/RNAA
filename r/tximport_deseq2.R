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

sanitize_name <- function(value) {
  gsub("[^A-Za-z0-9_.-]", "_", value)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  samplesheet_path <- required_arg(args, "samplesheet")
  quant_root <- required_arg(args, "quant-root")
  engine <- required_arg(args, "engine")
  tx2gene_path <- required_arg(args, "tx2gene")
  design_raw <- required_arg(args, "design")
  transform <- required_arg(args, "transform")
  counts_from_abundance <- required_arg(args, "counts-from-abundance")
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

  tx2gene <- utils::read.delim(tx2gene_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(tx2gene) < 2L) {
    stop("tx2gene must contain at least two columns", call. = FALSE)
  }
  tx2gene <- tx2gene[, 1:2]
  colnames(tx2gene) <- c("TXNAME", "GENEID")

  txi <- tximport::tximport(
    files,
    type = "kallisto",
    tx2gene = tx2gene,
    countsFromAbundance = counts_from_abundance
  )

  design_formula <- stats::as.formula(design_raw)
  dds <- DESeq2::DESeqDataSetFromTximport(
    txi = txi,
    colData = samplesheet,
    design = design_formula
  )
  dds <- DESeq2::DESeq(dds)

  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }

  gene_counts <- DESeq2::counts(dds, normalized = FALSE)
  gene_norm <- DESeq2::counts(dds, normalized = TRUE)
  utils::write.table(
    cbind(gene_id = rownames(gene_counts), as.data.frame(gene_counts)),
    file = file.path(outdir, "gene_counts.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    cbind(gene_id = rownames(gene_norm), as.data.frame(gene_norm)),
    file = file.path(outdir, "gene_norm_counts.tsv"),
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
  utils::write.table(
    cbind(gene_id = rownames(vst_mat), as.data.frame(vst_mat)),
    file = file.path(outdir, "vst.tsv"),
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
      gene_norm_counts = file.path(outdir, "gene_norm_counts.tsv"),
      vst = file.path(outdir, "vst.tsv"),
      vst_rds = file.path(outdir, "vst.rds"),
      size_factors = file.path(outdir, "size_factors.tsv"),
      de_tables = de_outputs
    ),
    session_info = utils::capture.output(sessionInfo())
  )
  dir.create(dirname(manifest_out), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(manifest, path = manifest_out, pretty = TRUE, auto_unbox = TRUE)
}

main()
