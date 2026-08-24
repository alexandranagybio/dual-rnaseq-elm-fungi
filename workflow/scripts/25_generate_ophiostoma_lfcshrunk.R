#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(readr)
  library(dplyr)
})

dds_path <- paste0(
  "results/ophiostoma/deseq2_results/objects/",
  "ophiostoma_dds_deseq2_fitted.rds"
)

raw_dir <- "results/ophiostoma/deseq2_results/tables"
out_dir <- raw_dir

alpha_threshold <- 0.05

if (!file.exists(dds_path)) {
  stop("Missing fitted Ophiostoma DESeq2 object: ", dds_path)
}

dds_self <- readRDS(dds_path)

if (!inherits(dds_self, "DESeqDataSet")) {
  stop("Input object is not a DESeqDataSet.")
}

expected_self_coefs <- c(
  "condition_interaction_vs_self",
  "condition_onu_vs_self"
)

missing_self_coefs <- setdiff(
  expected_self_coefs,
  resultsNames(dds_self)
)

if (length(missing_self_coefs) > 0) {
  stop(
    "Missing expected coefficients: ",
    paste(missing_self_coefs, collapse = ", ")
  )
}

read_raw_table <- function(contrast_name) {
  path <- file.path(
    raw_dir,
    paste0(contrast_name, "_all_genes.tsv")
  )

  if (!file.exists(path)) {
    stop("Missing canonical raw table: ", path)
  }

  read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE
  )
}

make_shrunk_table <- function(
    shrink_dds,
    coefficient,
    contrast_name
) {
  canonical_raw <- read_raw_table(contrast_name)

  named_raw <- results(
    shrink_dds,
    name = coefficient,
    alpha = alpha_threshold,
    independentFiltering = TRUE,
    cooksCutoff = TRUE
  )

  shrunk <- lfcShrink(
    shrink_dds,
    coef = coefficient,
    res = named_raw,
    type = "apeglm"
  )

  shrunk_df <- as.data.frame(shrunk) %>%
    tibble::rownames_to_column("gene_id") %>%
    transmute(
      gene_id,
      shrunk_log2FoldChange = log2FoldChange,
      shrunk_lfcSE = lfcSE
    )

  output <- canonical_raw %>%
    left_join(
      shrunk_df,
      by = "gene_id"
    )

  if (nrow(output) != nrow(canonical_raw)) {
    stop(
      contrast_name,
      ": row count changed after joining shrinkage estimates."
    )
  }

  if (anyDuplicated(output$gene_id) > 0) {
    stop(
      contrast_name,
      ": duplicated gene IDs after joining."
    )
  }

  if (any(is.na(output$shrunk_log2FoldChange))) {
    stop(
      contrast_name,
      ": missing shrunken LFC estimates."
    )
  }

  output
}

# ----------------------------------------------------------------------
# Direct coefficients in the canonical self-reference model
# ----------------------------------------------------------------------

interaction_vs_self <- make_shrunk_table(
  dds_self,
  "condition_interaction_vs_self",
  "interaction_vs_self"
)

onu_vs_self <- make_shrunk_table(
  dds_self,
  "condition_onu_vs_self",
  "onu_vs_self"
)

# ----------------------------------------------------------------------
# Reparameterize only for interaction versus ONU
#
# Existing size factors and dispersion estimates are retained.
# nbinomWaldTest recalculates coefficients under ONU as reference.
# ----------------------------------------------------------------------

dds_onu <- readRDS(dds_path)

dds_onu$condition <- relevel(
  factor(
    as.character(dds_onu$condition)
  ),
  ref = "onu"
)

design(dds_onu) <- ~ condition

dds_onu <- nbinomWaldTest(
  dds_onu,
  betaPrior = FALSE
)

required_onu_coef <- "condition_interaction_vs_onu"

if (!required_onu_coef %in% resultsNames(dds_onu)) {
  stop(
    "Missing reparameterized coefficient: ",
    required_onu_coef,
    "\nAvailable coefficients: ",
    paste(resultsNames(dds_onu), collapse = ", ")
  )
}

interaction_vs_onu <- make_shrunk_table(
  dds_onu,
  required_onu_coef,
  "interaction_vs_onu"
)

outputs <- list(
  interaction_vs_self = interaction_vs_self,
  interaction_vs_onu = interaction_vs_onu,
  onu_vs_self = onu_vs_self
)

summary_rows <- lapply(
  names(outputs),
  function(contrast_name) {
    x <- outputs[[contrast_name]]

    tibble(
      contrast = contrast_name,
      genes = nrow(x),
      significant_padj_0.05 = sum(
        !is.na(x$padj) &
          x$padj < 0.05
      ),
      raw_abs_lfc_gt_1 = sum(
        !is.na(x$padj) &
          x$padj < 0.05 &
          abs(x$log2FoldChange) > 1
      ),
      shrunk_abs_lfc_gt_1 = sum(
        !is.na(x$padj) &
          x$padj < 0.05 &
          abs(x$shrunk_log2FoldChange) > 1
      )
    )
  }
)

summary_table <- bind_rows(summary_rows)

for (contrast_name in names(outputs)) {
  write_tsv(
    outputs[[contrast_name]],
    file.path(
      out_dir,
      paste0(
        contrast_name,
        "_lfcshrunk.tsv"
      )
    )
  )
}

audit_dir <- "results/publication/deseq2_symmetry_audit"

dir.create(
  audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write_tsv(
  summary_table,
  file.path(
    audit_dir,
    "ophiostoma_lfc_shrinkage_summary.tsv"
  )
)

cat("\n============================================================\n")
cat("OPHIOSTOMA APEGLM SHRINKAGE COMPLETE\n")
cat("============================================================\n\n")

print(
  summary_table,
  n = Inf,
  width = Inf
)

cat("\nOutput files:\n")

for (contrast_name in names(outputs)) {
  cat(
    "  ",
    file.path(
      out_dir,
      paste0(
        contrast_name,
        "_lfcshrunk.tsv"
      )
    ),
    "\n",
    sep = ""
  )
}
