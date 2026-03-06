#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("missing required R package: jsonlite", call. = FALSE)
  }
}))

parse_args <- function(raw_args) {
  args <- list()
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
    args[[name]] <- value
    i <- i + 2L
  }
  args
}

required_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || identical(value, "")) {
    stop(paste("missing required argument --", name, sep = ""), call. = FALSE)
  }
  value
}

run_cmd_capture <- function(cmd, args) {
  out <- tryCatch(
    system2(cmd, args = args, stdout = TRUE, stderr = TRUE),
    error = function(e) character()
  )
  if (length(out) == 0) {
    return("")
  }
  paste(out, collapse = "\n")
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  run_id <- required_arg(args, "run")
  fastq1 <- required_arg(args, "fastq1")
  index <- required_arg(args, "index")
  outdir <- required_arg(args, "outdir")
  manifest_out <- required_arg(args, "manifest-out")
  threads_raw <- required_arg(args, "threads")
  reference_manifest <- args[["reference-manifest"]]
  fastq2 <- args[["fastq2"]]

  threads <- suppressWarnings(as.integer(threads_raw))
  if (is.na(threads) || threads < 1L) {
    stop("--threads must be a positive integer", call. = FALSE)
  }

  if (!file.exists(fastq1)) {
    stop(paste("fastq1 does not exist:", fastq1), call. = FALSE)
  }
  if (!is.null(fastq2) && !file.exists(fastq2)) {
    stop(paste("fastq2 does not exist:", fastq2), call. = FALSE)
  }
  if (!file.exists(index)) {
    stop(paste("kallisto index does not exist:", index), call. = FALSE)
  }
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }

  started_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  paired <- !is.null(fastq2) && !identical(fastq2, "")

  kallisto_args <- c(
    "quant",
    "-i", index,
    "-o", outdir,
    "-t", as.character(threads)
  )
  assumptions <- character()
  if (paired) {
    kallisto_args <- c(kallisto_args, fastq1, fastq2)
  } else {
    # Conservative single-end defaults until fragment params are user-configurable.
    kallisto_args <- c(kallisto_args, "--single", "-l", "200", "-s", "20", fastq1)
    assumptions <- c(assumptions, "single_end_defaults_l200_s20")
  }

  status <- system2("kallisto", args = kallisto_args, stdout = TRUE, stderr = TRUE)
  exit_status <- attr(status, "status")
  if (is.null(exit_status)) {
    exit_status <- 0L
  }
  if (exit_status != 0L) {
    msg <- paste(status, collapse = "\n")
    if (identical(msg, "")) {
      msg <- paste("kallisto exited with status", exit_status)
    }
    stop(msg, call. = FALSE)
  }

  abundance_h5 <- file.path(outdir, "abundance.h5")
  run_info_json <- file.path(outdir, "run_info.json")
  if (!file.exists(abundance_h5) || file.info(abundance_h5)$size <= 0) {
    stop("kallisto output missing abundance.h5", call. = FALSE)
  }
  if (!file.exists(run_info_json) || file.info(run_info_json)$size <= 0) {
    stop("kallisto output missing run_info.json", call. = FALSE)
  }

  dir.create(dirname(manifest_out), recursive = TRUE, showWarnings = FALSE)
  manifest <- list(
    schema_version = 1L,
    stage = "quant",
    id = run_id,
    started_at = started_at,
    finished_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    parameters = list(
      paired = paired,
      threads = threads,
      kallisto_args = kallisto_args,
      reference_manifest = reference_manifest
    ),
    tool_versions = list(
      R = R.version.string,
      kallisto = run_cmd_capture("kallisto", c("version"))
    ),
    assumptions = assumptions,
    session_info = utils::capture.output(sessionInfo())
  )
  jsonlite::write_json(manifest, path = manifest_out, pretty = TRUE, auto_unbox = TRUE)
}

main()
