#!/usr/bin/env Rscript
# Independent workflow oracle: execution of public DESeq2 APIs, never RNAA output.
stopifnot(getRversion() == "4.6.2", packageVersion("DESeq2") == "1.52.0")
options(digits = 17)
Sys.setlocale("LC_COLLATE", "C")
args <- commandArgs(TRUE)
out <- args[[1]]
input <- args[[2]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_matrix <- function(x, path) {
  write.table(data.frame(gene_id = rownames(x), x, check.names = FALSE),
              path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
original <- as.matrix(read.delim(input, row.names = 1, check.names = FALSE))
# Low expression, all-zero and a single extreme outlier exercise NA/filter policies.
counts <- rbind(original, all_zero = rep(0L, 8),
                low_expression = c(0, 1, 0, 0, 0, 0, 1, 0),
                extreme_outlier = c(1000000, 10, 12, 8, 9, 11, 10, 12))
for (case in c("categorical", "numeric", "numeric_interaction", "factor_interaction", "three_levels")) {
  directory <- file.path(out, case)
  dir.create(directory, showWarnings = FALSE)
  values <- counts
  if (case == "three_levels") {
    values <- cbind(values, counts[, 1:4])
    # A gene expressed only in C must have the character contrast A/B set to zero.
    values <- rbind(values, only_C = c(rep(0L, 8), 100L, 120L, 90L, 110L))
  }
  samples <- sprintf("sample%02d", seq_len(ncol(values)))
  colnames(values) <- samples
  metadata <- data.frame(run_accession = samples,
                         condition = factor(rep(c("A", "B", "C"), each = 4)[seq_along(samples)]),
                         age = rep(c(-1.5, -0.5, 0.5, 1.5), length.out = length(samples)),
                         batch = factor(rep(c("X", "X", "Y", "Y"), length.out = length(samples))),
                         row.names = samples)
  formula <- switch(case, categorical = ~ condition, numeric = ~ age + condition,
                    numeric_interaction = ~ age * condition,
                    factor_interaction = ~ batch * condition, three_levels = ~ condition)
  write_matrix(values, file.path(directory, "counts.tsv"))
  write.table(metadata, file.path(directory, "metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  dds <- DESeq2::DESeqDataSetFromMatrix(values, metadata, formula)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  for (direction in c("B_vs_A", "A_vs_B")) {
    levels <- strsplit(direction, "_vs_", fixed = TRUE)[[1]]
    results <- DESeq2::results(dds, contrast = c("condition", levels))
    write_matrix(as.data.frame(results), file.path(directory, paste0(direction, ".tsv")))
  }
  write_matrix(DESeq2::counts(dds, normalized = TRUE), file.path(directory, "normalized_counts.tsv"))
  write_matrix(SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)),
               file.path(directory, "vst.tsv"))
  write.table(data.frame(sample = samples, size_factor = DESeq2::sizeFactors(dds)),
              file.path(directory, "size_factors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(c(paste0("design\t", paste(deparse(formula), collapse = "")),
               paste0("fit_type\t", attr(DESeq2::dispersionFunction(dds), "fitType"))),
              file.path(directory, "fit.tsv"))
}
write.table(data.frame(package = c("R", "DESeq2"),
                        version = c(as.character(getRversion()), as.character(packageVersion("DESeq2")))),
            file.path(out, "versions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
