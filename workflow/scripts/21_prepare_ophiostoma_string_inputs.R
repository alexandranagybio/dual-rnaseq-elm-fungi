#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

# -------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------

deg_file <- file.path(
  "results",
  "ophiostoma",
  "deseq2_results",
  "tables",
  "interaction_vs_self_all_genes.tsv"
)

annotation_file <- file.path(
  "results",
  "ophiostoma",
  "annotation",
  "ophiostoma_gene_annotation_map.tsv"
)

protein_fasta <- file.path(
  "data",
  "reference",
  "ophiostoma",
  "ophiostoma.proteins.clean.fa"
)

functional_file <- file.path(
  "results",
  "ophiostoma",
  "functional_annotation",
  "ophiostoma_functional_annotation.tsv"
)

output_dir <- file.path(
  "results",
  "ophiostoma",
  "string_analysis"
)

input_dir <- file.path(output_dir, "input")
diagnostics_dir <- file.path(output_dir, "diagnostics")
tables_dir <- file.path(output_dir, "tables")

dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required file: ", path)
  }
}

read_fasta_ids <- function(path) {
  headers <- readLines(path, warn = FALSE)
  headers <- headers[str_starts(headers, ">")]

  cleaned_headers <- headers |>
    str_remove("^>")

  split_headers <- str_split_fixed(
    cleaned_headers,
    "\\s+",
    2
  )

  tibble(
    fasta_id = split_headers[, 1]
  )
}

normalise_missing <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "-", "\\N", "NA", "NaN")] <- NA_character_
  x
}

# -------------------------------------------------------------------------
# Validate inputs
# -------------------------------------------------------------------------

invisible(lapply(
  c(deg_file, annotation_file, protein_fasta, functional_file),
  stop_if_missing
))

# -------------------------------------------------------------------------
# Read inputs
# -------------------------------------------------------------------------

deg <- read_tsv(
  deg_file,
  show_col_types = FALSE,
  na = c("", "NA", "NaN")
)

annotation <- read_tsv(
  annotation_file,
  show_col_types = FALSE,
  na = c("", "NA", "NaN", "-", "\\N")
)

functional <- read_tsv(
  functional_file,
  show_col_types = FALSE,
  na = c("", "NA", "NaN", "-", "\\N")
)

fasta_ids <- read_fasta_ids(protein_fasta)

# -------------------------------------------------------------------------
# Structural validation
# -------------------------------------------------------------------------

required_deg_columns <- c(
  "gene_id",
  "contrast",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)

required_annotation_columns <- c(
  "gene_id",
  "protein_id",
  "transcript_id",
  "jgi_name",
  "portal_id",
  "go_count",
  "go_ids",
  "go_names",
  "kegg_count",
  "kegg_ecnum",
  "kegg_definition",
  "kegg_pathway",
  "kegg_pathway_class"
)

required_functional_columns <- c(
  "gene_id",
  "mrna_id",
  "signalp_is_sp",
  "dbcan_any_hit",
  "dbcan_high_confidence"
)

missing_deg_columns <- setdiff(required_deg_columns, names(deg))
missing_annotation_columns <- setdiff(
  required_annotation_columns,
  names(annotation)
)
missing_functional_columns <- setdiff(
  required_functional_columns,
  names(functional)
)

if (length(missing_deg_columns) > 0) {
  stop(
    "Missing DESeq2 columns: ",
    paste(missing_deg_columns, collapse = ", ")
  )
}

if (length(missing_annotation_columns) > 0) {
  stop(
    "Missing annotation columns: ",
    paste(missing_annotation_columns, collapse = ", ")
  )
}

if (length(missing_functional_columns) > 0) {
  stop(
    "Missing functional annotation columns: ",
    paste(missing_functional_columns, collapse = ", ")
  )
}

if (anyDuplicated(deg$gene_id)) {
  stop("Duplicated gene_id values found in DESeq2 table.")
}

if (anyDuplicated(annotation$gene_id)) {
  stop("Duplicated gene_id values found in annotation map.")
}

if (anyDuplicated(functional$gene_id)) {
  stop("Duplicated gene_id values found in functional annotation.")
}

if (anyDuplicated(fasta_ids$fasta_id)) {
  stop("Duplicated FASTA identifiers found.")
}

# -------------------------------------------------------------------------
# Harmonise identifiers
# -------------------------------------------------------------------------

annotation <- annotation |>
  mutate(
    protein_id = normalise_missing(protein_id),
    transcript_id = normalise_missing(transcript_id),
    jgi_name = normalise_missing(jgi_name),
    fasta_id = paste0("mRNA_", str_remove(gene_id, "^gene_"))
  )

functional <- functional |>
  mutate(
    mrna_id = normalise_missing(mrna_id)
  )

# Verify the explicit transcript/gene naming convention instead of assuming it.
gene_transcript_check <- annotation |>
  transmute(
    gene_id,
    protein_id,
    jgi_transcript_id = transcript_id,
    fasta_id,
    functional_mrna_id = functional$mrna_id[
      match(gene_id, functional$gene_id)
    ],
    fasta_id_present = fasta_id %in% fasta_ids$fasta_id,
    functional_matches_fasta =
      functional_mrna_id == fasta_id
  )

# -------------------------------------------------------------------------
# Join complete analysis table
# -------------------------------------------------------------------------

string_input <- deg |>
  left_join(
    annotation |>
      select(
        gene_id,
        protein_id,
        transcript_id,
        fasta_id,
        jgi_name,
        portal_id,
        go_count,
        go_ids,
        go_names,
        kegg_count,
        kegg_ecnum,
        kegg_definition,
        kegg_pathway,
        kegg_pathway_class
      ),
    by = "gene_id"
  ) |>
  left_join(
    functional |>
      select(
        gene_id,
        mrna_id,
        signalp_is_sp,
        dbcan_any_hit,
        dbcan_high_confidence,
        dbcan_n_tools,
        dbcan_recommended,
        dbcan_substrate
      ),
    by = "gene_id"
  ) |>
  mutate(
    is_tested = !is.na(padj),
    is_significant = is_tested & padj < 0.05,
    is_strong = is_significant & abs(log2FoldChange) > 1,
    direction = case_when(
      is_significant & log2FoldChange > 0 ~ "up",
      is_significant & log2FoldChange < 0 ~ "down",
      TRUE ~ "not_significant"
    ),
    strong_direction = case_when(
      is_strong & log2FoldChange > 0 ~ "up",
      is_strong & log2FoldChange < 0 ~ "down",
      TRUE ~ "not_strong"
    ),
    has_protein_id = !is.na(protein_id),
    has_transcript_id = !is.na(transcript_id),
    has_fasta_id = !is.na(fasta_id),
    fasta_id_in_fasta = fasta_id %in% fasta_ids$fasta_id,
    string_candidate_jgi = jgi_name,
    string_candidate_numeric = protein_id,
    string_candidate_transcript = transcript_id,
    string_candidate_fasta = fasta_id
  ) |>
  arrange(gene_id)

# -------------------------------------------------------------------------
# Write master table
# -------------------------------------------------------------------------

write_tsv(
  string_input,
  file.path(tables_dir, "ophiostoma_string_master_table.tsv"),
  na = ""
)

# -------------------------------------------------------------------------
# Write STRING candidate lists
# -------------------------------------------------------------------------

write_identifier_list <- function(data, column, path) {
  values <- data |>
    filter(!is.na(.data[[column]])) |>
    distinct(.data[[column]]) |>
    pull(.data[[column]])

  write_lines(values, path)
}

tested <- string_input |>
  filter(is_tested)

significant <- string_input |>
  filter(is_significant)

strong <- string_input |>
  filter(is_strong)

up <- string_input |>
  filter(is_significant, log2FoldChange > 0)

down <- string_input |>
  filter(is_significant, log2FoldChange < 0)

strong_up <- string_input |>
  filter(is_strong, log2FoldChange > 0)

strong_down <- string_input |>
  filter(is_strong, log2FoldChange < 0)

# Full tabular inputs

write_tsv(
  tested,
  file.path(input_dir, "interaction_vs_self_tested_background.tsv"),
  na = ""
)

write_tsv(
  significant,
  file.path(input_dir, "interaction_vs_self_significant.tsv"),
  na = ""
)

write_tsv(
  strong,
  file.path(input_dir, "interaction_vs_self_strong.tsv"),
  na = ""
)

write_tsv(
  up,
  file.path(input_dir, "interaction_vs_self_significant_up.tsv"),
  na = ""
)

write_tsv(
  down,
  file.path(input_dir, "interaction_vs_self_significant_down.tsv"),
  na = ""
)

write_tsv(
  strong_up,
  file.path(input_dir, "interaction_vs_self_strong_up.tsv"),
  na = ""
)

write_tsv(
  strong_down,
  file.path(input_dir, "interaction_vs_self_strong_down.tsv"),
  na = ""
)

# Plain identifier files for testing STRING recognition.

for (dataset_name in c(
  "tested",
  "significant",
  "strong",
  "up",
  "down",
  "strong_up",
  "strong_down"
)) {
  dataset <- get(dataset_name)

  write_identifier_list(
    dataset,
    "jgi_name",
    file.path(
      input_dir,
      paste0(dataset_name, "_jgi_identifiers.txt")
    )
  )

  write_identifier_list(
    dataset,
    "protein_id",
    file.path(
      input_dir,
      paste0(dataset_name, "_numeric_protein_ids.txt")
    )
  )

  write_identifier_list(
    dataset,
    "transcript_id",
    file.path(
      input_dir,
      paste0(dataset_name, "_transcript_ids.txt")
    )
  )
}

# -------------------------------------------------------------------------
# Diagnostics
# -------------------------------------------------------------------------

mapping_summary <- tibble(
  metric = c(
    "deseq2_rows",
    "tested_genes",
    "significant_genes_padj_lt_0.05",
    "strong_genes_padj_lt_0.05_abs_lfc_gt_1",
    "significant_up",
    "significant_down",
    "strong_up",
    "strong_down",
    "annotation_rows",
    "protein_fasta_entries",
    "tested_with_protein_id",
    "tested_with_transcript_id",
    "tested_fasta_id_present_in_fasta",
    "significant_with_protein_id",
    "strong_with_protein_id",
    "tested_with_go",
    "tested_with_kegg",
    "significant_with_go",
    "significant_with_kegg"
  ),
  value = c(
    nrow(deg),
    nrow(tested),
    nrow(significant),
    nrow(strong),
    nrow(up),
    nrow(down),
    nrow(strong_up),
    nrow(strong_down),
    nrow(annotation),
    nrow(fasta_ids),
    sum(tested$has_protein_id),
    sum(tested$has_transcript_id),
    sum(tested$fasta_id_in_fasta),
    sum(significant$has_protein_id),
    sum(strong$has_protein_id),
    sum(!is.na(tested$go_ids)),
    sum(!is.na(tested$kegg_pathway)),
    sum(!is.na(significant$go_ids)),
    sum(!is.na(significant$kegg_pathway))
  )
) |>
  mutate(
    expected = c(
      8560,
      8560,
      NA,
      NA,
      NA,
      NA,
      NA,
      NA,
      8640,
      8640,
      8560,
      8560,
      8560,
      NA,
      NA,
      NA,
      NA,
      NA,
      NA
    ),
    status = case_when(
      is.na(expected) ~ "INFO",
      value == expected ~ "PASS",
      TRUE ~ "FAIL"
    )
  )

write_tsv(
  mapping_summary,
  file.path(diagnostics_dir, "identifier_mapping_summary.tsv"),
  na = ""
)

write_tsv(
  gene_transcript_check,
  file.path(diagnostics_dir, "gene_transcript_fasta_validation.tsv"),
  na = ""
)

write_tsv(
  string_input |>
    filter(
      is_tested &
        (
          is.na(protein_id) |
            is.na(fasta_id) |
            !fasta_id_in_fasta
        )
    ),
  file.path(diagnostics_dir, "tested_genes_with_mapping_problem.tsv"),
  na = ""
)

write_tsv(
  anti_join(
    fasta_ids,
    annotation |> select(fasta_id),
    by = "fasta_id"
  ),
  file.path(diagnostics_dir, "fasta_ids_without_annotation.tsv"),
  na = ""
)

write_tsv(
  anti_join(
    annotation |> select(fasta_id),
    fasta_ids,
    by = "fasta_id"
  ),
  file.path(diagnostics_dir, "annotation_ids_without_fasta.tsv"),
  na = ""
)

# -------------------------------------------------------------------------
# Console report
# -------------------------------------------------------------------------

cat("\nOphiostoma STRING input preparation completed.\n\n")

print(mapping_summary, n = Inf)

cat("\nOutput directory:\n  ", output_dir, "\n", sep = "")

if (any(mapping_summary$status == "FAIL")) {
  stop(
    "\nOne or more structural validation checks failed. ",
    "Inspect the diagnostics before proceeding to STRING."
  )
}

cat("\nPASS: core mapping validation completed successfully.\n")
