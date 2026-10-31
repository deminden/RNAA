#!/usr/bin/env Rscript
# Independent expected values: public tximport and DESeq2 APIs only.
stopifnot(getRversion() == "4.6.2", packageVersion("tximport") == "1.40.0",
          packageVersion("DESeq2") == "1.52.0")
options(digits = 17)
Sys.setlocale("LC_COLLATE", "C")
out <- commandArgs(TRUE)[[1]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_matrix <- function(x, path) {
  write.table(data.frame(gene_id = rownames(x), x, check.names = FALSE),
              path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
# RNAA's declared GENCODE identity policy, independent R implementation.
canonical_id <- function(x) {
  x <- sub("\\|.*$", "", x)
  suffix <- ifelse(grepl("_PAR_Y$", x), "_PAR_Y", "")
  paste0(sub("\\.[0-9]+$", "", sub("_PAR_Y$", "", x)), suffix)
}
ids <- c("ENST1.2|ENSG1.5|iso1|", "ENST2.3|ENSG1.5|iso2|",
         "ENST3.2_PAR_Y|ENSG2.7_PAR_Y|label|", "ENST4.4|ENSG3.8|label|",
         "ENST5.1|ENSG4.2|zero|", "UNMAPPED.1|no_gene|")
genes <- c("ENSG1.5", "ENSG1.5", "ENSG2.7_PAR_Y", "ENSG3.8", "ENSG4.2")
count <- rbind(c(10.25, 0, 7.5, 1.5), c(5.25, 30.5, 0, 8.5),
               c(100.5, 80.5, 90.5, 120.5), c(0, 6.5, 4.5, 12.5),
               c(0, 0, 0, 0), c(20, 20, 20, 20))
length <- rbind(c(100, 120, 140, 160), c(500, 550, 600, 650),
                c(200, 220, 240, 260), c(100, 150, 250, 350),
                c(100, 200, 300, 400), c(100, 100, 100, 100))
samples <- paste0("sample", 1:4)
mapping <- data.frame(transcript = ids[1:5], gene = genes)
for (case in c("fractional_and_zero_abundance", "zero_length_without_mass")) {
  directory <- file.path(out, case)
  dir.create(directory, showWarnings = FALSE)
  effective <- length
  if (case == "zero_length_without_mass") effective[4, 1] <- 0
  files <- setNames(file.path(directory, paste0(samples, ".tsv")), samples)
  write.table(mapping, file.path(directory, "tx2gene.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  for (sample in seq_along(samples)) {
    rate <- ifelse(count[, sample] == 0, 0, count[, sample] / effective[, sample])
    tpm <- rate / sum(rate) * 1e6
    write.table(data.frame(target_id = ids, length = rep(1000, length(ids)),
                           eff_length = effective[, sample], est_counts = count[, sample], tpm = tpm),
                files[sample], sep = "\t", quote = FALSE, row.names = FALSE)
  }
  reader <- function(file) {
    data <- read.delim(file, check.names = FALSE)
    data$target_id <- canonical_id(data$target_id)
    data
  }
  map <- data.frame(transcript = canonical_id(mapping$transcript), gene = canonical_id(mapping$gene))
  for (mode in c("no", "scaledTPM", "lengthScaledTPM")) {
    imported <- tximport::tximport(files, type = "kallisto", tx2gene = map,
                                  countsFromAbundance = mode, importer = reader, dropInfReps = TRUE)
    write_matrix(round(imported$counts), file.path(directory, paste0(mode, "_counts.tsv")))
    if (mode == "no") {
      length_factors <- imported$length / exp(rowMeans(log(imported$length)))
      write_matrix(length_factors, file.path(directory, "length_factors.tsv"))
      dds <- DESeq2::DESeqDataSetFromTximport(imported,
              data.frame(row.names = samples, dummy = seq_along(samples)), ~ 1)
      dds <- DESeq2::estimateSizeFactors(dds)
      write_matrix(DESeq2::normalizationFactors(dds), file.path(directory, "normalization_factors.tsv"))
      write_matrix(DESeq2::counts(dds, normalized = TRUE), file.path(directory, "normalized_counts.tsv"))
    }
  }
}
write.table(data.frame(package = c("R", "tximport", "DESeq2"),
                        version = c(as.character(getRversion()), as.character(packageVersion("tximport")),
                                    as.character(packageVersion("DESeq2")))),
            file.path(out, "versions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
