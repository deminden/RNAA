#!/usr/bin/env Rscript
stopifnot(getRversion() == "4.6.2", packageVersion("WGCNA") == "1.74",
          packageVersion("HellCor") == "1.3")
options(digits = 17)
Sys.setlocale("LC_COLLATE", "C")
set.seed(20261030)
out <- commandArgs(TRUE)[[1]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_matrix <- function(x, path) {
  write.table(data.frame(gene_id = rownames(x), x, check.names = FALSE), path,
              sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
cases <- list(
  ties_and_outliers = rbind(ties = c(1, 1, 2, 3, 4, 4, 5, 7),
       inverse = 10:3, outlier = c(1, 2, 1, 2, 1, 2, 1, 50),
       zero_mad = c(1, 1, 1, 1, 1, 1, 1, 9), constant = rep(3, 8),
       nonlinear = c(-4, -1, 0, 1, 4, 9, 16, 25)),
  minimum_two_samples = rbind(up = c(1, 2), down = c(3, 1), constant = c(4, 4),
                             other = c(9, 0)),
  three_samples = rbind(up = c(1, 2, 3), tied = c(1, 1, 4), down = c(5, 3, 0),
                       constant = c(2, 2, 2)),
  near_perfect = rbind(up = 1:8, almost_up = c(1:7, 8.00001), inverse = 8:1,
                       irregular = c(1, 3, 4, 2, 6, 5, 8, 7), constant = rep(2, 8))
)
for (name in names(cases)) {
  directory <- file.path(out, name)
  dir.create(directory, showWarnings = FALSE)
  data <- cases[[name]]
  colnames(data) <- paste0("sample", seq_len(ncol(data)))
  write_matrix(data, file.path(directory, "input.tsv"))
  constant <- apply(data, 1, function(row) length(unique(row)) == 1)
  for (method in c("pearson", "spearman", "kendall", "bicor", "hellcor")) {
    result <- if (method == "bicor") {
      suppressWarnings(WGCNA::bicor(t(data), use = "all.obs", maxPOutliers = 1,
                                    pearsonFallback = "individual", nThreads = 1))
    } else if (method == "hellcor") {
      matrix <- matrix(NA_real_, nrow(data), nrow(data), dimnames = list(rownames(data), rownames(data)))
      for (i in seq_len(nrow(data))) for (j in i:nrow(data)) {
        value <- suppressWarnings(HellCor::HellCor(data[i, ], data[j, ], Kmax = 20L,
                       Lmax = 20L, alpha = 6, pval.comp = FALSE, C.version = TRUE)$Hcor)
        matrix[i, j] <- matrix[j, i] <- value
      }
      matrix
    } else suppressWarnings(stats::cor(t(data), method = method))
    write_matrix(result, file.path(directory, paste0(method, "_raw.tsv")))
    # Matrix contract: undefined constant-row pairs stay missing, while variable
    # self-correlations are one. Raw oracle estimates remain separately auditable.
    result[constant, ] <- NA_real_
    result[, constant] <- NA_real_
    diag(result)[!constant] <- 1
    endpoints <- is.finite(result) & abs(abs(result) - 1) <= 8 * .Machine$double.eps
    result[endpoints] <- sign(result[endpoints])
    write_matrix(result, file.path(directory, paste0(method, ".tsv")))
    # Derive topology only from R coefficients, independently of native writers.
    write_edges <- function(edges, filename) {
      rows <- if (length(edges)) do.call(rbind, edges) else data.frame(gene_a = character(), gene_b = character(), r = numeric())
      write.table(rows, file.path(directory, paste0(method, "_", filename, ".tsv")),
                  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    }
    edges <- list()
    for (i in seq_len(nrow(result))) {
      targets <- which(is.finite(result[i, ]) & seq_len(nrow(result)) != i)
      targets <- head(targets[order(-abs(result[i, targets]), targets)], 2)
      for (j in targets) edges[[length(edges) + 1]] <- data.frame(gene_a = rownames(result)[i], gene_b = rownames(result)[j], r = result[i, j])
    }
    write_edges(edges, "topk2")
    for (threshold in c(0, 0.75, 1)) {
      edges <- list()
      for (i in seq_len(nrow(result))) for (j in seq_len(nrow(result))) {
        if (j > i && is.finite(result[i, j]) && abs(result[i, j]) >= threshold)
          edges[[length(edges) + 1]] <- data.frame(gene_a = rownames(result)[i], gene_b = rownames(result)[j], r = result[i, j])
      }
      write_edges(edges, paste0("threshold", threshold))
    }
  }
}
write.table(data.frame(package = c("R", "WGCNA", "HellCor"),
               version = c(as.character(getRversion()), as.character(packageVersion("WGCNA")),
                           as.character(packageVersion("HellCor")))),
            file.path(out, "versions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
