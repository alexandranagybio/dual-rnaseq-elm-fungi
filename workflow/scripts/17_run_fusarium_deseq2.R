#!/usr/bin/env Rscript

# ---------------------------------------------------------------------------
# 17_run_fusarium_deseq2.R
#
# Purpose:
#   Import Fusarium Salmon transcript estimates with tximport, aggregate
#   Trinity transcripts to Trinity genes, and run gene-level DESeq2 analysis.
#
# Canonical inputs:
#   - config/samples.tsv
#   - data/external/fusarium_salmon/quants/
#   - data/external/fusarium_assembly/
#       Fusarium_pure.Trinity.fasta.gene_trans_map
#
# Comparison:
#   Fusarium interaction versus Fusarium self-control
#
# Filtering:
#   Retain genes with counts >= 10 in at least 2 samples.
#
# Usage:
#   Rscript workflow/scripts/17_run_fusarium_deseq2.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(apeglm)
  library(readr)
})

# ---------------------------------------------------------------------------
# Resolve repository paths
# ---------------------------------------------------------------------------

command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command_args, value = TRUE)

if (length(file_arg) != 1L) {
  stop("Could not determine the path of this script.")
}

script_path <- normalizePath(
  sub("^--file=", "", file_arg),
  mustWork = TRUE
)

repo_root <- normalizePath(
  file.path(dirname(script_path), "../.."),
  mustWork = TRUE
)

sample_tsv <- file.path(
  repo_root,
  "config",
  "samples.tsv"
)

quant_root <- file.path(
  repo_root,
  "data",
  "external",
  "fusarium_salmon",
  "quants"
)

gene_trans_map <- file.path(
  repo_root,
  "data",
  "external",
  "fusarium_assembly",
  "Fusarium_pure.Trinity.fasta.gene_trans_map"
)

out_dir <- file.path(
  repo_root,
  "results",
  "fusarium",
  "deseq2_results"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ---------------------------------------------------------------------------
# Validate canonical inputs
# ---------------------------------------------------------------------------

required_paths <- c(
  sample_tsv,
  quant_root,
  gene_trans_map
)

missing_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_paths) > 0L) {
  stop(
    "Missing required Fusarium input path(s):\n",
    paste(missing_paths, collapse = "\n")
  )
}

# ---------------------------------------------------------------------------
# Read Fusarium sample metadata
# ---------------------------------------------------------------------------

samples_all <- read_tsv(
  sample_tsv,
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
)

required_columns <- c(
  "sample_id",
  "condition",
  "biological_replicate",
  "include_fusarium"
)

missing_columns <- setdiff(
  required_columns,
  colnames(samples_all)
)

if (length(missing_columns) > 0L) {
  stop(
    "Missing required column(s) from config/samples.tsv: ",
    paste(missing_columns, collapse = ", ")
  )
}

samples <- samples_all[
  samples_all$include_fusarium == "yes",
  required_columns
]

samples$analysis_condition <- ifelse(
  samples$condition == "fusarium_interaction",
  "interaction",
  ifelse(
    samples$condition == "fusarium_self",
    "control",
    NA_character_
  )
)

if (anyNA(samples$analysis_condition)) {
  stop(
    "Unexpected Fusarium condition label(s): ",
    paste(
      unique(samples$condition[is.na(samples$analysis_condition)]),
      collapse = ", "
    )
  )
}

expected_samples <- c(
  "146",
  "147",
  "148",
  "158",
  "159",
  "160"
)

if (!identical(sort(samples$sample_id), expected_samples)) {
  stop(
    "Expected Fusarium samples ",
    paste(expected_samples, collapse = ", "),
    " but found ",
    paste(sort(samples$sample_id), collapse = ", ")
  )
}

if (anyDuplicated(samples$sample_id)) {
  stop("Duplicated Fusarium sample IDs were found.")
}

samples$analysis_condition <- factor(
  samples$analysis_condition,
  levels = c("control", "interaction")
)

# ---------------------------------------------------------------------------
# Resolve Salmon quantification files
# ---------------------------------------------------------------------------

files <- file.path(
  quant_root,
  paste0("P21526_", samples$sample_id),
  "quant.sf"
)

names(files) <- samples$sample_id

missing_files <- files[!file.exists(files)]

if (length(missing_files) > 0L) {
  stop(
    "Missing Salmon quantification file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

empty_files <- files[file.info(files)$size == 0]

if (length(empty_files) > 0L) {
  stop(
    "Empty Salmon quantification file(s):\n",
    paste(empty_files, collapse = "\n")
  )
}

# ---------------------------------------------------------------------------
# Read and validate Trinity gene-transcript map
#
# Trinity gene_trans_map orientation:
#   column 1 = gene
#   column 2 = transcript
#
# tximport requires:
#   column 1 = transcript
#   column 2 = gene
# ---------------------------------------------------------------------------

gene_trans_raw <- read.delim(
  gene_trans_map,
  header = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(gene_trans_raw) < 2L) {
  stop("Trinity gene_trans_map contains fewer than two columns.")
}

gene_trans_raw <- gene_trans_raw[, 1:2]

colnames(gene_trans_raw) <- c(
  "gene",
  "transcript"
)

if (anyNA(gene_trans_raw) ||
    any(gene_trans_raw$gene == "") ||
    any(gene_trans_raw$transcript == "")) {
  stop("Missing identifiers found in Trinity gene_trans_map.")
}

if (!all(grepl("_i[0-9]+$", gene_trans_raw$transcript))) {
  stop(
    "The second gene_trans_map column does not consistently contain ",
    "Trinity transcript IDs ending in _i<number>."
  )
}

if (anyDuplicated(gene_trans_raw$transcript)) {
  stop("Duplicated transcript IDs found in Trinity gene_trans_map.")
}

tx2gene <- gene_trans_raw[, c("transcript", "gene")]

message(
  "Validated Trinity map: ",
  format(nrow(tx2gene), big.mark = ","),
  " transcripts; ",
  format(length(unique(tx2gene$gene)), big.mark = ","),
  " genes."
)

# ---------------------------------------------------------------------------
# Import Salmon estimates and construct DESeq2 dataset
# ---------------------------------------------------------------------------

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene
)

col_data <- data.frame(
  condition = samples$analysis_condition,
  biological_replicate = samples$biological_replicate,
  row.names = samples$sample_id,
  stringsAsFactors = FALSE
)

if (!identical(colnames(txi$counts), rownames(col_data))) {
  stop(
    "Sample order differs between tximport counts and sample metadata."
  )
}

dds <- DESeqDataSetFromTximport(
  txi,
  colData = col_data,
  design = ~ condition
)

input_gene_count <- nrow(dds)

keep <- rowSums(counts(dds) >= 10) >= 2

dds <- dds[keep, ]

retained_gene_count <- nrow(dds)
filtered_gene_count <- input_gene_count - retained_gene_count

if (retained_gene_count == 0L) {
  stop("No genes remained after the low-count filter.")
}

message(
  "DESeq2 genes before filtering: ",
  format(input_gene_count, big.mark = ",")
)

message(
  "DESeq2 genes retained: ",
  format(retained_gene_count, big.mark = ",")
)

message(
  "DESeq2 genes removed: ",
  format(filtered_gene_count, big.mark = ",")
)

# ---------------------------------------------------------------------------
# Run DESeq2
# ---------------------------------------------------------------------------

dds <- DESeq(dds)

coef_name <- "condition_interaction_vs_control"

if (!coef_name %in% resultsNames(dds)) {
  stop(
    "Expected DESeq2 coefficient not found: ",
    coef_name,
    "\nAvailable coefficients: ",
    paste(resultsNames(dds), collapse = ", ")
  )
}

res_raw <- results(
  dds,
  contrast = c(
    "condition",
    "interaction",
    "control"
  ),
  alpha = 0.05
)

res_shrunk <- lfcShrink(
  dds,
  coef = coef_name,
  type = "apeglm"
)

raw_table <- cbind(
  gene_id = rownames(as.data.frame(res_raw)),
  as.data.frame(res_raw)
)

shrunk_table <- cbind(
  gene_id = rownames(as.data.frame(res_shrunk)),
  as.data.frame(res_shrunk)
)

write.table(
  raw_table,
  file = file.path(
    out_dir,
    "fusarium_interaction_vs_self_raw.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  shrunk_table,
  file = file.path(
    out_dir,
    "fusarium_interaction_vs_self_lfcshrunk.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  dds,
  file.path(
    out_dir,
    "fusarium_deseq2_dataset.rds"
  )
)

# ---------------------------------------------------------------------------
# QC plots
#
# blind = TRUE is used for unsupervised sample-level QC.
# ---------------------------------------------------------------------------

vsd <- vst(
  dds,
  blind = TRUE
)

pdf(
  file.path(
    out_dir,
    "fusarium_pca_vst_blind.pdf"
  )
)

print(
  plotPCA(
    vsd,
    intgroup = "condition"
  )
)

dev.off()

pdf(
  file.path(
    out_dir,
    "fusarium_ma_lfcshrunk.pdf"
  )
)

plotMA(
  res_shrunk,
  ylim = c(-6, 6)
)

dev.off()

# ---------------------------------------------------------------------------
# Summary and provenance
# ---------------------------------------------------------------------------

significant <- sum(
  res_shrunk$padj < 0.05,
  na.rm = TRUE
)

upregulated <- sum(
  res_shrunk$padj < 0.05 &
    res_shrunk$log2FoldChange > 0,
  na.rm = TRUE
)

downregulated <- sum(
  res_shrunk$padj < 0.05 &
    res_shrunk$log2FoldChange < 0,
  na.rm = TRUE
)

summary_lines <- c(
  paste("repository_root", ".", sep = "\t"),
  paste(
    "sample_table",
    "config/samples.tsv",
    sep = "\t"
  ),
  paste(
    "quantification_root",
    "data/external/fusarium_salmon/quants",
    sep = "\t"
  ),
  paste(
    "gene_transcript_map",
    paste0(
      "data/external/fusarium_assembly/",
      "Fusarium_pure.Trinity.fasta.gene_trans_map"
    ),
    sep = "\t"
  ),
  paste("sample_count", nrow(samples), sep = "\t"),
  paste(
    "sample_ids",
    paste(samples$sample_id, collapse = ","),
    sep = "\t"
  ),
  paste(
    "reference_condition",
    "control",
    sep = "\t"
  ),
  paste(
    "contrast",
    "interaction_vs_control",
    sep = "\t"
  ),
  paste(
    "input_gene_count",
    input_gene_count,
    sep = "\t"
  ),
  paste(
    "retained_gene_count",
    retained_gene_count,
    sep = "\t"
  ),
  paste(
    "filtered_gene_count",
    filtered_gene_count,
    sep = "\t"
  ),
  paste(
    "filter",
    "counts >= 10 in at least 2 samples",
    sep = "\t"
  ),
  paste(
    "significant_padj_lt_0.05",
    significant,
    sep = "\t"
  ),
  paste(
    "up_in_interaction",
    upregulated,
    sep = "\t"
  ),
  paste(
    "down_in_interaction",
    downregulated,
    sep = "\t"
  ),
  paste(
    "vst_blind",
    "TRUE",
    sep = "\t"
  )
)

writeLines(
  summary_lines,
  file.path(
    out_dir,
    "fusarium_deseq2_run_summary.tsv"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(
    out_dir,
    "fusarium_deseq2_session_info.txt"
  )
)

message(
  "Significant genes (padj < 0.05): ",
  significant
)

message(
  "Upregulated in interaction: ",
  upregulated
)

message(
  "Downregulated in interaction: ",
  downregulated
)

message(
  "SUCCESS: Fusarium DESeq2 analysis completed."
)

message(
  "Output directory: ",
  out_dir
)
