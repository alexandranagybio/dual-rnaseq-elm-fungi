#!/usr/bin/env Rscript

# 06_build_ophiostoma_deseq2_dataset.R
#
# Build and validate the Ophiostoma gene-level DESeq2 dataset without yet
# performing inferential comparisons. This separates input validation from
# biological analysis.
#
# Usage:
#   Rscript workflow/scripts/06_build_ophiostoma_deseq2_dataset.R
#
# Required R package:
#   DESeq2

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) != 1L) {
  stop("Could not determine the script path.", call. = FALSE)
}

script_path <- normalizePath(sub("^--file=", "", file_arg))
script_dir <- dirname(script_path)
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

input_dir <- file.path(repo_root, "results", "ophiostoma", "gene_counts")
matrix_path <- file.path(input_dir, "ophiostoma_gene_counts_matrix.tsv")
metadata_path <- file.path(input_dir, "sample_metadata.tsv")
output_rds <- file.path(input_dir, "ophiostoma_gene_level_dds_unfiltered.rds")
session_path <- file.path(input_dir, "06_build_ophiostoma_deseq2_dataset.sessionInfo.txt")
summary_path <- file.path(input_dir, "06_build_ophiostoma_deseq2_dataset.validation.tsv")

expected_genes <- 8640L
expected_samples <- as.character(149:157)

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("R package 'DESeq2' is not installed in this environment.", call. = FALSE)
}

if (!file.exists(matrix_path) || file.info(matrix_path)$size == 0) {
  stop("Missing or empty count matrix: ", matrix_path, call. = FALSE)
}

if (!file.exists(metadata_path) || file.info(metadata_path)$size == 0) {
  stop("Missing or empty metadata: ", metadata_path, call. = FALSE)
}

counts_df <- read.delim(
  matrix_path,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

metadata <- read.delim(
  metadata_path,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (nrow(counts_df) != expected_genes) {
  stop(
    "Expected ", expected_genes, " genes, observed ", nrow(counts_df), ".",
    call. = FALSE
  )
}

if (anyDuplicated(rownames(counts_df))) {
  stop("Duplicated gene IDs in count matrix.", call. = FALSE)
}

if (!identical(colnames(counts_df), expected_samples)) {
  stop(
    "Count-matrix columns are not in canonical order 149-157.",
    call. = FALSE
  )
}

if (!identical(rownames(metadata), expected_samples)) {
  stop(
    "Metadata rows are not in canonical order 149-157.",
    call. = FALSE
  )
}

if (!identical(colnames(counts_df), rownames(metadata))) {
  stop("Count columns and metadata rows do not match exactly.", call. = FALSE)
}

count_matrix <- as.matrix(counts_df)

if (anyNA(count_matrix)) {
  stop("NA values detected in count matrix.", call. = FALSE)
}

if (any(count_matrix < 0)) {
  stop("Negative values detected in count matrix.", call. = FALSE)
}

if (any(count_matrix != round(count_matrix))) {
  stop("Non-integer values detected in count matrix.", call. = FALSE)
}

storage.mode(count_matrix) <- "integer"

metadata$condition <- factor(
  metadata$condition,
  levels = c("self", "onu", "interaction")
)
metadata$replicate <- factor(metadata$replicate)

if (anyNA(metadata$condition)) {
  stop("Unexpected condition label in metadata.", call. = FALSE)
}

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ condition
)

saveRDS(dds, output_rds)

validation <- data.frame(
  metric = c(
    "gene_rows",
    "samples",
    "zero_total_genes",
    "design",
    "condition_reference",
    "dataset_path"
  ),
  value = c(
    nrow(dds),
    ncol(dds),
    sum(rowSums(DESeq2::counts(dds)) == 0L),
    "~ condition",
    levels(metadata$condition)[1],
    output_rds
  ),
  stringsAsFactors = FALSE
)

write.table(
  validation,
  file = summary_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

sink(session_path)
sessionInfo()
sink()

cat("PASS: DESeq2 dataset created and validated.\n")
cat("Genes:", format(nrow(dds), big.mark = ","), "\n")
cat("Samples:", ncol(dds), "\n")
cat("Design: ~ condition\n")
cat("Reference condition:", levels(metadata$condition)[1], "\n")
cat("Dataset:", output_rds, "\n")
cat("Validation:", summary_path, "\n")
cat("Session info:", session_path, "\n")
