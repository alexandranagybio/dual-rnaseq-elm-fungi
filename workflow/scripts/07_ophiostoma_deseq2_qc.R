#!/usr/bin/env Rscript

# 07_ophiostoma_deseq2_qc.R
#
# Purpose:
#   Perform pre-inference quality control on the validated Ophiostoma
#   gene-level DESeq2 dataset.
#
# This script:
#   1. Removes genes with zero total counts.
#   2. Estimates DESeq2 size factors.
#   3. Applies variance-stabilizing transformation (VST).
#   4. Produces library-size, size-factor, PCA, sample-distance,
#      and per-sample correlation outputs.
#   5. Saves the filtered but not yet statistically fitted DESeq2 object.
#
# It does NOT:
#   - run DESeq();
#   - calculate differential-expression results;
#   - choose an adjusted-p-value or fold-change threshold;
#   - remove any sample automatically.
#
# Usage:
#   Rscript workflow/scripts/07_ophiostoma_deseq2_qc.R

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) != 1L) {
  stop("Could not determine the script path.", call. = FALSE)
}

script_path <- normalizePath(sub("^--file=", "", file_arg))
script_dir <- dirname(script_path)
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

input_dir <- file.path(repo_root, "results", "ophiostoma", "gene_counts")
output_dir <- file.path(repo_root, "results", "ophiostoma", "deseq2_qc")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dds_path <- file.path(
  input_dir,
  "ophiostoma_gene_level_dds_unfiltered.rds"
)

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("R package 'DESeq2' is not installed.", call. = FALSE)
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("R package 'ggplot2' is not installed.", call. = FALSE)
}

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  stop("R package 'pheatmap' is not installed.", call. = FALSE)
}

if (!file.exists(dds_path)) {
  stop("Missing DESeq2 object: ", dds_path, call. = FALSE)
}

dds <- readRDS(dds_path)

expected_samples <- as.character(149:157)

if (nrow(dds) != 8640L) {
  stop("Expected 8,640 input genes, observed ", nrow(dds), ".", call. = FALSE)
}

if (!identical(colnames(dds), expected_samples)) {
  stop("Unexpected sample order in DESeq2 object.", call. = FALSE)
}

raw_counts <- DESeq2::counts(dds)
keep_nonzero <- rowSums(raw_counts) > 0L
dds_filtered <- dds[keep_nonzero, ]

if (sum(!keep_nonzero) != 80L) {
  warning(
    "Expected 80 all-zero genes from the prior validation, observed ",
    sum(!keep_nonzero),
    "."
  )
}

dds_filtered <- DESeq2::estimateSizeFactors(dds_filtered)
size_factors <- DESeq2::sizeFactors(dds_filtered)

if (anyNA(size_factors) || any(!is.finite(size_factors)) || any(size_factors <= 0)) {
  stop("Invalid DESeq2 size factors were estimated.", call. = FALSE)
}

vst <- DESeq2::vst(dds_filtered, blind = TRUE)
vst_matrix <- SummarizedExperiment::assay(vst)

sample_metadata <- as.data.frame(SummarizedExperiment::colData(dds_filtered))
sample_metadata$sample_id <- rownames(sample_metadata)

library_sizes <- colSums(raw_counts)

qc_table <- data.frame(
  sample_id = expected_samples,
  condition = as.character(sample_metadata[expected_samples, "condition"]),
  replicate = as.character(sample_metadata[expected_samples, "replicate"]),
  raw_library_size = as.numeric(library_sizes[expected_samples]),
  size_factor = as.numeric(size_factors[expected_samples]),
  normalized_library_size = as.numeric(
    colSums(DESeq2::counts(dds_filtered, normalized = TRUE))[expected_samples]
  ),
  stringsAsFactors = FALSE
)

write.table(
  qc_table,
  file = file.path(output_dir, "sample_qc_metrics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    metric = c(
      "input_genes",
      "all_zero_genes_removed",
      "genes_retained_for_qc",
      "samples",
      "vst_blind"
    ),
    value = c(
      nrow(dds),
      sum(!keep_nonzero),
      nrow(dds_filtered),
      ncol(dds_filtered),
      TRUE
    )
  ),
  file = file.path(output_dir, "qc_validation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -------------------------------------------------------------------------
# Library size plot
# -------------------------------------------------------------------------

library_plot <- ggplot2::ggplot(
  qc_table,
  ggplot2::aes(
    x = sample_id,
    y = raw_library_size,
    fill = condition
  )
) +
  ggplot2::geom_col() +
  ggplot2::scale_y_continuous(labels = scales::label_comma()) +
  ggplot2::labs(
    title = "Assigned Ophiostoma gene-count library sizes",
    x = "Sample",
    y = "Assigned fragments"
  ) +
  ggplot2::theme_bw(base_size = 12)

ggplot2::ggsave(
  filename = file.path(output_dir, "library_sizes.pdf"),
  plot = library_plot,
  width = 7,
  height = 5
)

# -------------------------------------------------------------------------
# Size factor plot
# -------------------------------------------------------------------------

size_factor_plot <- ggplot2::ggplot(
  qc_table,
  ggplot2::aes(
    x = sample_id,
    y = size_factor,
    fill = condition
  )
) +
  ggplot2::geom_col() +
  ggplot2::geom_hline(yintercept = 1, linetype = 2) +
  ggplot2::labs(
    title = "DESeq2 size factors",
    x = "Sample",
    y = "Size factor"
  ) +
  ggplot2::theme_bw(base_size = 12)

ggplot2::ggsave(
  filename = file.path(output_dir, "size_factors.pdf"),
  plot = size_factor_plot,
  width = 7,
  height = 5
)

# -------------------------------------------------------------------------
# PCA
# -------------------------------------------------------------------------

pca_data <- DESeq2::plotPCA(
  vst,
  intgroup = "condition",
  returnData = TRUE
)

percent_variance <- round(
  100 * attr(pca_data, "percentVar"),
  digits = 1
)

pca_data$sample_id <- rownames(pca_data)

write.table(
  pca_data,
  file = file.path(output_dir, "pca_coordinates.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pca_plot <- ggplot2::ggplot(
  pca_data,
  ggplot2::aes(
    x = PC1,
    y = PC2,
    colour = condition,
    label = sample_id
  )
) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_text(
    nudge_y = 0.5,
    check_overlap = TRUE
  ) +
  ggplot2::xlab(paste0("PC1: ", percent_variance[1], "% variance")) +
  ggplot2::ylab(paste0("PC2: ", percent_variance[2], "% variance")) +
  ggplot2::ggtitle("Ophiostoma sample PCA: blind VST") +
  ggplot2::theme_bw(base_size = 12)

ggplot2::ggsave(
  filename = file.path(output_dir, "pca_blind_vst.pdf"),
  plot = pca_plot,
  width = 7,
  height = 5.5
)

# -------------------------------------------------------------------------
# Sample distances and correlations
# -------------------------------------------------------------------------

sample_dist <- dist(t(vst_matrix))
sample_dist_matrix <- as.matrix(sample_dist)

write.table(
  sample_dist_matrix,
  file = file.path(output_dir, "sample_distance_matrix.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

annotation_col <- data.frame(
  condition = sample_metadata[colnames(vst_matrix), "condition"],
  row.names = colnames(vst_matrix)
)

pdf(
  file.path(output_dir, "sample_distance_heatmap.pdf"),
  width = 7,
  height = 6.5
)
pheatmap::pheatmap(
  sample_dist_matrix,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  main = "Sample distances: blind VST"
)
dev.off()

correlation_matrix <- stats::cor(vst_matrix, method = "pearson")

write.table(
  correlation_matrix,
  file = file.path(output_dir, "sample_pearson_correlations.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

pdf(
  file.path(output_dir, "sample_correlation_heatmap.pdf"),
  width = 7,
  height = 6.5
)
pheatmap::pheatmap(
  correlation_matrix,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  main = "Sample correlations: blind VST"
)
dev.off()

# Save objects needed for reproducible follow-up analysis.
saveRDS(
  dds_filtered,
  file.path(output_dir, "ophiostoma_gene_level_dds_nonzero.rds")
)

saveRDS(
  vst,
  file.path(output_dir, "ophiostoma_gene_level_vst_blind.rds")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "07_ophiostoma_deseq2_qc.sessionInfo.txt")
)

cat("PASS: pre-inference DESeq2 QC completed.\n")
cat("Input genes:", nrow(dds), "\n")
cat("All-zero genes removed:", sum(!keep_nonzero), "\n")
cat("Genes retained:", nrow(dds_filtered), "\n")
cat("Samples:", ncol(dds_filtered), "\n")
cat("Output directory:", output_dir, "\n")
cat("\nSample QC metrics:\n")
print(qc_table, row.names = FALSE)
