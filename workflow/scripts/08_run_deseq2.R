#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 120)
suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
input_rds <- if (length(args) >= 1L) args[[1]] else "results/ophiostoma/deseq2_qc/ophiostoma_dds_filtered.rds"
output_dir <- if (length(args) >= 2L) args[[2]] else "results/ophiostoma/deseq2_results"

expected_gene_count <- 8560L
expected_sample_count <- 9L
expected_conditions <- c("self", "onu", "interaction")
alpha_threshold <- 0.05
strong_lfc_threshold <- 1

for (d in c(output_dir,
            file.path(output_dir, "tables"),
            file.path(output_dir, "diagnostics"),
            file.path(output_dir, "objects"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

tables_dir <- file.path(output_dir, "tables")
diagnostics_dir <- file.path(output_dir, "diagnostics")
objects_dir <- file.path(output_dir, "objects")
log_file <- file.path(output_dir, "08_run_deseq2.log")

log_message <- function(...) {
  msg <- paste0(...)
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

fail <- function(...) {
  msg <- paste0(...)
  log_message("ERROR: ", msg)
  stop(msg, call. = FALSE)
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = TRUE, na = "NA")
}

sanitize_result_table <- function(res, contrast_name) {
  df <- as.data.frame(res)
  df$gene_id <- rownames(df)
  rownames(df) <- NULL
  df$contrast <- contrast_name
  df[, c("gene_id", "contrast", "baseMean", "log2FoldChange", "lfcSE",
         "stat", "pvalue", "padj")]
}

summarize_contrast <- function(df, contrast_name, numerator, denominator) {
  significant <- !is.na(df$padj) & df$padj < alpha_threshold
  strong <- significant & !is.na(df$log2FoldChange) &
    abs(df$log2FoldChange) > strong_lfc_threshold

  data.frame(
    contrast = contrast_name,
    numerator = numerator,
    denominator = denominator,
    total_rows = nrow(df),
    genes_tested = sum(!is.na(df$pvalue)),
    genes_with_padj = sum(!is.na(df$padj)),
    genes_with_na_padj = sum(is.na(df$padj)),
    significant_padj_lt_0_05 = sum(significant),
    significant_up = sum(significant & df$log2FoldChange > 0, na.rm = TRUE),
    significant_down = sum(significant & df$log2FoldChange < 0, na.rm = TRUE),
    strong_padj_lt_0_05_abs_lfc_gt_1 = sum(strong),
    strong_up = sum(strong & df$log2FoldChange > strong_lfc_threshold, na.rm = TRUE),
    strong_down = sum(strong & df$log2FoldChange < -strong_lfc_threshold, na.rm = TRUE),
    alpha = alpha_threshold,
    strong_abs_log2fc_threshold = strong_lfc_threshold,
    stringsAsFactors = FALSE
  )
}

if (!file.exists(input_rds)) fail("Input DESeqDataSet not found: ", input_rds)
log_message("Reading filtered DESeqDataSet: ", input_rds)
dds <- readRDS(input_rds)

if (!inherits(dds, "DESeqDataSet")) fail("Input object is not a DESeqDataSet.")
if (nrow(dds) != expected_gene_count) fail("Expected ", expected_gene_count, " genes; observed ", nrow(dds), ".")
if (ncol(dds) != expected_sample_count) fail("Expected ", expected_sample_count, " samples; observed ", ncol(dds), ".")
if (!"condition" %in% colnames(colData(dds))) fail("colData(dds) lacks a condition column.")

observed_conditions <- unique(as.character(colData(dds)$condition))
if (!setequal(observed_conditions, expected_conditions)) {
  fail("Expected conditions: ", paste(expected_conditions, collapse = ", "),
       "; observed: ", paste(observed_conditions, collapse = ", "), ".")
}

dds$condition <- factor(as.character(dds$condition), levels = c("self", "onu", "interaction"))
design(dds) <- ~ condition

if (levels(dds$condition)[1] != "self") fail("Reference level is not self.")
if (any(rowSums(counts(dds)) == 0L)) fail("Filtered object still contains all-zero genes.")
if (anyDuplicated(rownames(dds)) > 0L) fail("Duplicate gene identifiers detected.")
if (anyDuplicated(colnames(dds)) > 0L) fail("Duplicate sample identifiers detected.")
if (!all(colnames(dds) == rownames(colData(dds)))) fail("Count columns and metadata rows are not identically ordered.")

write_tsv(data.frame(
  sample = colnames(dds),
  condition = as.character(dds$condition),
  raw_library_size = colSums(counts(dds)),
  stringsAsFactors = FALSE
), file.path(tables_dir, "sample_overview.tsv"))

log_message("Validated input: ", nrow(dds), " genes, ", ncol(dds),
            " samples, design = ~ condition, reference = self.")

set.seed(1)
log_message("Running DESeq().")
dds <- DESeq(dds, test = "Wald", fitType = "parametric", sfType = "ratio",
             minReplicatesForReplace = 7, quiet = FALSE)
log_message("DESeq() completed.")

non_converged <- if ("betaConv" %in% colnames(mcols(dds))) {
  sum(!mcols(dds)$betaConv, na.rm = TRUE)
} else NA_integer_

if (!is.na(non_converged) && non_converged > 0L) {
  log_message("WARNING: ", non_converged, " genes have betaConv = FALSE.")
}

if (any(!is.finite(dispersions(dds))) || any(dispersions(dds) <= 0, na.rm = TRUE)) {
  fail("Non-finite or non-positive dispersion estimates detected.")
}
if (any(!is.finite(sizeFactors(dds))) || any(sizeFactors(dds) <= 0)) {
  fail("Non-finite or non-positive size factors detected.")
}

fitted_dds_path <- file.path(objects_dir, "ophiostoma_dds_deseq2_fitted.rds")
saveRDS(dds, fitted_dds_path)

write_tsv(data.frame(
  metric = c("input_genes", "samples", "conditions", "reference_level", "design",
             "genes_with_nonconverged_beta", "minimum_dispersion", "median_dispersion",
             "maximum_dispersion", "minimum_size_factor", "median_size_factor", "maximum_size_factor"),
  value = c(nrow(dds), ncol(dds), paste(levels(dds$condition), collapse = ","),
            levels(dds$condition)[1], "~ condition", non_converged,
            min(dispersions(dds)), median(dispersions(dds)), max(dispersions(dds)),
            min(sizeFactors(dds)), median(sizeFactors(dds)), max(sizeFactors(dds))),
  stringsAsFactors = FALSE
), file.path(diagnostics_dir, "model_diagnostics.tsv"))

write_tsv(data.frame(
  gene_id = rownames(dds),
  baseMean = mcols(dds)$baseMean,
  dispersion = dispersions(dds),
  maxCooks = if ("maxCooks" %in% colnames(mcols(dds))) mcols(dds)$maxCooks else NA_real_,
  betaConverged = if ("betaConv" %in% colnames(mcols(dds))) mcols(dds)$betaConv else NA,
  stringsAsFactors = FALSE
), file.path(diagnostics_dir, "gene_model_diagnostics.tsv"))

contrasts <- list(
  interaction_vs_self = c("condition", "interaction", "self"),
  onu_vs_self = c("condition", "onu", "self"),
  interaction_vs_onu = c("condition", "interaction", "onu")
)

contrast_summaries <- list()
all_complete_results <- list()

for (contrast_name in names(contrasts)) {
  contrast_vector <- contrasts[[contrast_name]]
  numerator <- contrast_vector[[2]]
  denominator <- contrast_vector[[3]]
  log_message("Extracting ", contrast_name, ": ", numerator, " vs ", denominator, ".")

  res <- results(dds, contrast = contrast_vector, alpha = alpha_threshold,
                 independentFiltering = TRUE, cooksCutoff = TRUE)
  result_df <- sanitize_result_table(res, contrast_name)
  result_df <- result_df[match(rownames(dds), result_df$gene_id), ]

  significant_df <- result_df[!is.na(result_df$padj) & result_df$padj < alpha_threshold, , drop = FALSE]
  strong_df <- significant_df[!is.na(significant_df$log2FoldChange) &
                                abs(significant_df$log2FoldChange) > strong_lfc_threshold, , drop = FALSE]

  significant_df <- significant_df[order(significant_df$padj,
                                         -abs(significant_df$log2FoldChange),
                                         na.last = TRUE), , drop = FALSE]
  strong_df <- strong_df[order(strong_df$padj,
                               -abs(strong_df$log2FoldChange),
                               na.last = TRUE), , drop = FALSE]

  write_tsv(result_df, file.path(tables_dir, paste0(contrast_name, "_all_genes.tsv")))
  write_tsv(significant_df, file.path(tables_dir, paste0(contrast_name, "_significant_padj_0.05.tsv")))
  write_tsv(strong_df, file.path(tables_dir, paste0(contrast_name, "_strong_padj_0.05_abs_log2fc_gt_1.tsv")))

  summary_df <- summarize_contrast(result_df, contrast_name, numerator, denominator)
  contrast_summaries[[contrast_name]] <- summary_df
  all_complete_results[[contrast_name]] <- result_df

  filter_threshold <- if (!is.null(metadata(res)$filterThreshold)) metadata(res)$filterThreshold else NA
  write_tsv(data.frame(
    field = c("contrast", "numerator", "denominator", "alpha", "independent_filtering",
              "cooks_cutoff", "filter_threshold", "strong_abs_log2fc_threshold"),
    value = c(contrast_name, numerator, denominator, alpha_threshold, TRUE, TRUE,
              filter_threshold, strong_lfc_threshold),
    stringsAsFactors = FALSE
  ), file.path(diagnostics_dir, paste0(contrast_name, "_results_metadata.tsv")))

  log_message(contrast_name, ": significant = ", summary_df$significant_padj_lt_0_05,
              "; strong = ", summary_df$strong_padj_lt_0_05_abs_lfc_gt_1,
              "; NA padj = ", summary_df$genes_with_na_padj, ".")
}

summary_table <- do.call(rbind, contrast_summaries)
rownames(summary_table) <- NULL
write_tsv(summary_table, file.path(tables_dir, "contrast_summary.tsv"))

combined_results <- do.call(rbind, all_complete_results)
rownames(combined_results) <- NULL
write_tsv(combined_results, file.path(tables_dir, "all_contrasts_complete_results.tsv"))

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

write_tsv(data.frame(
  field = c("script", "run_timestamp", "input_rds", "output_dir", "R_version",
            "DESeq2_version", "design", "reference_level", "alpha",
            "strong_abs_log2fc_threshold", "fitType", "sfType", "test"),
  value = c("08_run_deseq2.R", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
            normalizePath(input_rds), normalizePath(output_dir), R.version.string,
            as.character(packageVersion("DESeq2")), "~ condition", "self",
            alpha_threshold, strong_lfc_threshold, "parametric", "ratio", "Wald"),
  stringsAsFactors = FALSE
), file.path(output_dir, "run_info.tsv"))

required_outputs <- c(
  fitted_dds_path,
  file.path(tables_dir, "contrast_summary.tsv"),
  file.path(tables_dir, "all_contrasts_complete_results.tsv"),
  file.path(diagnostics_dir, "model_diagnostics.tsv"),
  file.path(output_dir, "sessionInfo.txt"),
  file.path(output_dir, "run_info.tsv")
)
missing_outputs <- required_outputs[!file.exists(required_outputs)]
if (length(missing_outputs) > 0L) fail("Required outputs missing: ", paste(missing_outputs, collapse = ", "))

final_status <- if (!is.na(non_converged) && non_converged > 0L) "PASS_WITH_WARNINGS" else "PASS"

write_tsv(data.frame(
  check = c("input_is_DESeqDataSet", "gene_count_matches_expected", "sample_count_matches_expected",
            "condition_levels_match_expected", "reference_level_is_self", "no_all_zero_genes",
            "no_duplicate_gene_ids", "sample_order_matches_metadata", "positive_finite_size_factors",
            "positive_finite_dispersions", "required_outputs_created", "overall_status"),
  value = c(TRUE, nrow(dds) == expected_gene_count, ncol(dds) == expected_sample_count,
            setequal(levels(dds$condition), expected_conditions), levels(dds$condition)[1] == "self",
            !any(rowSums(counts(dds)) == 0L), anyDuplicated(rownames(dds)) == 0L,
            all(colnames(dds) == rownames(colData(dds))),
            all(is.finite(sizeFactors(dds))) && all(sizeFactors(dds) > 0),
            all(is.finite(dispersions(dds))) && all(dispersions(dds) > 0),
            length(missing_outputs) == 0L, final_status),
  stringsAsFactors = FALSE
), file.path(output_dir, "validation.tsv"))

cat("\n", final_status, "\n\n", sep = "")
cat("Genes fitted:", nrow(dds), "\n")
cat("Samples:", ncol(dds), "\n")
cat("Design: ~ condition\n")
cat("Reference: self\n")
cat("Non-converged gene coefficients:", ifelse(is.na(non_converged), "not available", non_converged), "\n\n")

for (i in seq_len(nrow(summary_table))) {
  cat(summary_table$contrast[[i]], ":\n", sep = "")
  cat("  Significant (padj < 0.05): ", summary_table$significant_padj_lt_0_05[[i]], "\n", sep = "")
  cat("  Strong (padj < 0.05 and |log2FC| > 1): ",
      summary_table$strong_padj_lt_0_05_abs_lfc_gt_1[[i]], "\n", sep = "")
  cat("  NA adjusted p-values: ", summary_table$genes_with_na_padj[[i]], "\n\n", sep = "")
}

cat("Results written to:", output_dir, "\n")
