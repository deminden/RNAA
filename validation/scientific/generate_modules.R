#!/usr/bin/env Rscript
# Independent PCA, graph and enrichment oracles; no native outputs are inputs.
stopifnot(getRversion() == "4.6.2", packageVersion("WGCNA") == "1.74",
          packageVersion("fgsea") == "1.38.0")
options(digits = 17)
Sys.setlocale("LC_COLLATE", "C")
out <- commandArgs(TRUE)[[1]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write_matrix <- function(x, path) write_table(data.frame(gene_id = rownames(x), x, check.names = FALSE), path)
score <- function(data) {
  standardized <- t(scale(t(data)))
  gram <- crossprod(standardized)
  eig <- eigen(gram, symmetric = TRUE)
  leading <- which(eig$values[1] - eig$values <= 64 * .Machine$double.eps * ncol(data) * eig$values[1])
  if (length(leading) == 1) {
    # Use the independent SVD PCA API for an ordinary leading component.
    pc <- prcomp(t(data), center = TRUE, scale. = TRUE)$x[, 1]
  } else {
    basis <- eig$vectors[, leading, drop = FALSE]
    for (axis in seq_len(ncol(data))) {
      pc <- as.numeric(basis %*% basis[axis, ])
      if (sqrt(sum(pc^2)) > 64 * .Machine$double.eps * ncol(data)) break
    }
  }
  pc <- pc / sqrt(sum(pc^2))
  average <- colSums(standardized)
  alignment <- sum(pc * average)
  if (abs(alignment) > 64 * .Machine$double.eps * sqrt(sum(average^2)) * sqrt(sum(pc^2))) {
    if (alignment < 0) pc <- -pc
  } else {
    largest <- max(abs(pc))
    anchor <- which(abs(pc) >= largest - largest * 64 * .Machine$double.eps * ncol(data))[[1]]
    if (pc[anchor] < 0) pc <- -pc
  }
  list(values = as.numeric(scale(pc)), variance_explained = eig$values[1] / sum(diag(gram)), multiplicity = length(leading))
}
cases <- list(
 unequal_scales = rbind(c(1, 2, 4, 3, 6, 5, 8, 9) * 1000,
                        c(3, 1, 4, 2, 5, 8, 6, 7) / 100,
                        c(0, 1, 0, 3, 1, 2, 4, 2),
                        c(4, 3, 5, 1, 2, 0, -1, -3)),
 zero_average = rbind(c(-3, -1, 1, 3), c(3, 1, -1, -3)),
 tied_leading = rbind(c(-1, 1, -1, 1), c(-1, -1, 1, 1)),
 three_samples = rbind(c(1, 4, 2), c(2, 1, 3), c(1, 2, 5))
)
for (name in names(cases)) {
  directory <- file.path(out, name); dir.create(directory, showWarnings = FALSE)
  data <- cases[[name]]; rownames(data) <- paste0("gene", seq_len(nrow(data))); colnames(data) <- paste0("sample", seq_len(ncol(data)))
  expected <- score(data)
  write_matrix(data, file.path(directory, "input.tsv"))
  write_table(data.frame(sample = colnames(data), score = expected$values), file.path(directory, "scores.tsv"))
  write_table(data.frame(variance_explained = expected$variance_explained, leading_multiplicity = expected$multiplicity), file.path(directory, "diagnostics.tsv"))
  if (expected$multiplicity == 1 && name != "zero_average") {
    raw <- WGCNA::moduleEigengenes(t(data), rep("blue", nrow(data)), scale = TRUE, impute = FALSE,
                                  subHubs = FALSE, trapErrors = FALSE)$eigengenes[[1]]
    write_table(data.frame(sample = colnames(data), raw_wgcna = raw, normalized_wgcna = as.numeric(scale(raw))), file.path(directory, "wgcna.tsv"))
  }
}
directory <- file.path(out, "workflow"); dir.create(directory, showWarnings = FALSE)
data <- sapply(1:8, function(sample) sapply(0:29, function(gene) {
  base <- if (gene < 15) sample else 9 - sample
  scale <- if (gene == 0) 1000 else if (gene == 1) 0.01 else (gene + 1) / 3
  (base + ((gene * sample + sample^2) %% 7 - 3) * 0.6) * scale
}))
rownames(data) <- sprintf("G%02d", 0:29); colnames(data) <- paste0("sample", 1:8)
data <- rbind(data, constant = rep(3, 8))
write_matrix(data, file.path(directory, "input.tsv"))
pathways <- list(A = sprintf("G%02d", 0:14), B = sprintf("G%02d", 15:29), MIX = sprintf("G%02d", c(0:4, 15:19)))
writeLines(c(paste(c("A", "synthetic", pathways$A, "UNMAPPED"), collapse = "\t"),
             paste(c("B", "synthetic", pathways$B), collapse = "\t"),
             paste(c("MIX", "synthetic", pathways$MIX), collapse = "\t"), "NO_MATCH\tsynthetic\tUNKNOWN"), file.path(directory, "pathways.gmt"))
corr <- suppressWarnings(cor(t(data), method = "spearman"))
corr[is.finite(corr) & abs(abs(corr) - 1) <= 8 * .Machine$double.eps] <- sign(corr[is.finite(corr) & abs(abs(corr) - 1) <= 8 * .Machine$double.eps])
neighbors <- lapply(seq_len(nrow(data)), function(i) {
  targets <- which(is.finite(corr[i, ]) & corr[i, ] > 0 & seq_len(nrow(data)) != i)
  head(targets[order(-corr[i, targets], targets)], min(20, nrow(data) - 1))
})
visited <- rep(FALSE, nrow(data)); modules <- list()
for (i in seq_len(nrow(data))) {
  if (visited[i]) next
  queue <- i; visited[i] <- TRUE; members <- integer()
  while (length(queue)) {
    node <- queue[[1]]; queue <- queue[-1]; members <- c(members, node)
    mutual <- neighbors[[node]][vapply(neighbors[[node]], function(j) node %in% neighbors[[j]], logical(1))]
    fresh <- mutual[!visited[mutual]]; visited[fresh] <- TRUE; queue <- c(queue, fresh)
  }
  if (length(members) >= 2) modules[[length(modules) + 1]] <- sort(members)
}
modules <- modules[order(-lengths(modules), vapply(modules, min, integer(1)))]
membership <- list(); scores <- list(); rankings <- list(); enrichment <- list()
for (index in seq_along(modules)) {
  id <- sprintf("M%04d", index); members <- modules[[index]]; pc <- score(data[members, , drop = FALSE])
  connectivity <- vapply(members, function(gene) sum(pmax(corr[gene, setdiff(members, gene)], 0)), numeric(1))
  ordered <- order(-round(connectivity, 10), members)
  membership[[index]] <- data.frame(module_id = id, module_size = length(members), gene = rownames(data)[members[ordered]], rank_in_module = seq_along(members), connectivity = connectivity[ordered])
  scores[[index]] <- data.frame(module_id = id, t(pc$values), check.names = FALSE); colnames(scores[[index]]) <- c("module_id", colnames(data))
  eligible <- which(apply(data, 1, function(row) length(unique(row)) > 1))
  r <- vapply(eligible, function(gene) cor(data[gene, ], pc$values, method = "spearman"), numeric(1))
  r[abs(abs(r) - 1) <= 8 * .Machine$double.eps] <- sign(r[abs(abs(r) - 1) <= 8 * .Machine$double.eps])
  capped <- pmax(-0.9999999, pmin(0.9999999, r)); t <- capped * sqrt((ncol(data) - 2) / (1 - capped^2))
  rankings[[index]] <- data.frame(module_id = id, gene = rownames(data)[eligible], spearman_r = r, t_value = t)
  set.seed(42)
  result <- suppressWarnings(fgsea::fgseaSimple(pathways, setNames(t, rownames(data)[eligible]), nperm = 100,
                                               minSize = 10, maxSize = length(eligible) - 1, nproc = 1))
  enrichment[[index]] <- data.frame(module_id = id, pathway = result$pathway, es = result$ES, nes = result$NES,
                                    p_value = result$pval, padj = result$padj, pathway_size = result$size,
                                    leading_edge_size = lengths(result$leadingEdge))
}
write_table(do.call(rbind, membership), file.path(directory, "membership.tsv"))
write_table(do.call(rbind, scores), file.path(directory, "scores.tsv"))
write_table(do.call(rbind, rankings), file.path(directory, "rankings.tsv"))
write_table(do.call(rbind, enrichment), file.path(directory, "enrichment.tsv"))
write_table(data.frame(package = c("R", "WGCNA", "fgsea"), version = c(as.character(getRversion()), as.character(packageVersion("WGCNA")), as.character(packageVersion("fgsea")))), file.path(out, "versions.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
