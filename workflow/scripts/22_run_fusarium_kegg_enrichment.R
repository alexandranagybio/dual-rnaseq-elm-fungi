#!/usr/bin/env Rscript

# 22_run_fusarium_kegg_enrichment.R
#
# Directional KEGG over-representation analysis for Fusarium salinense
# Trinity gene-level DESeq2 results.
#
# Primary analysis: padj < 0.05, split into up/down genes.
# Sensitivity analysis: padj < 0.05 and abs(log2FoldChange) > 1.
# Background: all 10,414 genes tested by DESeq2.
#
# KEGG annotation is supplied by:
#   results/fusarium/kegg_annotation/fusarium_kegg_term2gene.tsv
#
# Optional pathway names are read from:
#   results/fusarium/kegg_annotation/fusarium_kegg_pathway_metadata.tsv
#
# The current metadata may contain pathway IDs as names. Descriptive names can
# be added later without changing the enrichment calculations.

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ggplot2)
})

EXPECTED_TESTED_GENES <- 10414L
EXPECTED_STRUCTURAL_GENES <- 15192L
EXPECTED_KEGG_GENES <- 2503L
EXPECTED_GENE_PATHWAY_PAIRS <- 10897L
EXPECTED_PATHWAYS <- 383L

PADJ_CUTOFF <- 0.05
STRONG_LFC_CUTOFF <- 1
MIN_GENE_SET_SIZE <- 5L
MAX_GENE_SET_SIZE <- 500L

annotation_dir <- file.path("results", "fusarium", "annotation")
deseq2_dir <- file.path("results", "fusarium", "deseq2_results")
kegg_annotation_dir <- file.path(
  "results", "fusarium", "kegg_annotation"
)
output_dir <- file.path(
  "results", "fusarium", "kegg_enrichment"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gene_annotation_file <- file.path(
  file.path("results", "fusarium", "functional_annotation", "tables"), "fusarium_gene_functional_annotation.tsv"
)
term2gene_file <- file.path(
  kegg_annotation_dir, "fusarium_kegg_term2gene.tsv"
)
metadata_file <- file.path(
  kegg_annotation_dir, "fusarium_kegg_pathway_metadata.tsv"
)

contrast_files <- c(
  interaction_vs_self = file.path(
    deseq2_dir, "fusarium_interaction_vs_self_raw.tsv"
  )
)

stop_with <- function(...) stop(paste0(...), call. = FALSE)

read_tsv_checked <- function(path, required_columns = NULL) {
  if (!file.exists(path)) stop_with("Missing input file: ", path)

  x <- read.delim(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    na.strings = c("", "NA", "NaN"),
    stringsAsFactors = FALSE
  )

  if (!is.null(required_columns)) {
    missing_columns <- setdiff(required_columns, names(x))
    if (length(missing_columns) > 0L) {
      stop_with(
        "Input ", path, " is missing required columns: ",
        paste(missing_columns, collapse = ", ")
      )
    }
  }

  x
}

write_empty_tsv <- function(path) {
  write.table(
    data.frame(),
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

parse_ratio <- function(x) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)
  vapply(
    parts,
    function(z) {
      if (length(z) != 2L) return(NA_real_)
      a <- suppressWarnings(as.numeric(z[1]))
      b <- suppressWarnings(as.numeric(z[2]))
      if (is.na(a) || is.na(b) || b == 0) return(NA_real_)
      a / b
    },
    numeric(1)
  )
}

# Read and validate the structural gene table.
gene_annotation <- read_tsv_checked(
  gene_annotation_file,
  "gene_id"
)
gene_annotation$gene_id <- trimws(
  as.character(gene_annotation$gene_id)
)

if (any(is.na(gene_annotation$gene_id) | gene_annotation$gene_id == "")) {
  stop_with("The Fusarium gene annotation contains missing gene_id values.")
}
if (nrow(gene_annotation) != EXPECTED_STRUCTURAL_GENES) {
  stop_with(
    "Expected ", EXPECTED_STRUCTURAL_GENES,
    " structural genes, observed ", nrow(gene_annotation), "."
  )
}
if (anyDuplicated(gene_annotation$gene_id)) {
  stop_with("Duplicated gene_id values in Fusarium gene annotation.")
}

structural_genes <- sort(unique(gene_annotation$gene_id))

# Read and validate TERM2GENE.
term2gene_raw <- read_tsv_checked(
  term2gene_file,
  c("pathway", "gene_id")
)
term2gene_raw$pathway <- trimws(as.character(term2gene_raw$pathway))
term2gene_raw$gene_id <- trimws(as.character(term2gene_raw$gene_id))

invalid_pairs <- (
  is.na(term2gene_raw$pathway) |
    term2gene_raw$pathway == "" |
    is.na(term2gene_raw$gene_id) |
    term2gene_raw$gene_id == ""
)
if (any(invalid_pairs)) {
  stop_with("TERM2GENE contains missing pathway or gene_id values.")
}
if (anyDuplicated(term2gene_raw[, c("pathway", "gene_id")])) {
  stop_with("TERM2GENE contains duplicated pathway-gene pairs.")
}

invalid_pathway_ids <- !grepl(
  "^map[0-9]{5}$",
  term2gene_raw$pathway
)
if (any(invalid_pathway_ids)) {
  stop_with(
    "TERM2GENE contains invalid canonical KEGG pathway IDs: ",
    paste(
      head(unique(term2gene_raw$pathway[invalid_pathway_ids]), 20L),
      collapse = ", "
    )
  )
}

unmapped <- setdiff(
  unique(term2gene_raw$gene_id),
  structural_genes
)
if (length(unmapped) > 0L) {
  stop_with(
    "TERM2GENE genes absent from the structural annotation: ",
    paste(head(unmapped, 20L), collapse = ", ")
  )
}

if (length(unique(term2gene_raw$gene_id)) != EXPECTED_KEGG_GENES) {
  stop_with(
    "Expected ", EXPECTED_KEGG_GENES,
    " KEGG-annotated genes, observed ",
    length(unique(term2gene_raw$gene_id)), "."
  )
}
if (nrow(term2gene_raw) != EXPECTED_GENE_PATHWAY_PAIRS) {
  stop_with(
    "Expected ", EXPECTED_GENE_PATHWAY_PAIRS,
    " gene-pathway pairs, observed ", nrow(term2gene_raw), "."
  )
}
if (length(unique(term2gene_raw$pathway)) != EXPECTED_PATHWAYS) {
  stop_with(
    "Expected ", EXPECTED_PATHWAYS,
    " pathways, observed ",
    length(unique(term2gene_raw$pathway)), "."
  )
}

term2gene <- term2gene_raw
names(term2gene) <- c("term", "gene")

# Read optional pathway metadata and construct TERM2NAME.
metadata_available <- FALSE
term2name_used <- FALSE
metadata_unique <- NULL
term2name <- NULL

if (file.exists(metadata_file)) {
  pathway_metadata <- read_tsv_checked(
    metadata_file,
    c("pathway", "pathway_name")
  )
  pathway_metadata$pathway <- trimws(
    as.character(pathway_metadata$pathway)
  )
  pathway_metadata$pathway_name <- trimws(
    as.character(pathway_metadata$pathway_name)
  )

  if (any(is.na(pathway_metadata$pathway) |
          pathway_metadata$pathway == "")) {
    stop_with("KEGG pathway metadata contains missing pathway IDs.")
  }
  if (anyDuplicated(pathway_metadata$pathway)) {
    stop_with("KEGG pathway metadata contains duplicated pathway IDs.")
  }
  if (!setequal(
    unique(term2gene_raw$pathway),
    pathway_metadata$pathway
  )) {
    stop_with("Pathway sets differ between TERM2GENE and metadata.")
  }

  metadata_unique <- pathway_metadata
  valid_names <- (
    !is.na(metadata_unique$pathway_name) &
      metadata_unique$pathway_name != ""
  )

  if (all(valid_names)) {
    term2name <- metadata_unique[
      , c("pathway", "pathway_name"), drop = FALSE
    ]
    names(term2name) <- c("term", "name")
    term2name_used <- TRUE
  }

  metadata_available <- TRUE
}

run_enrichment <- function(
  selected_genes,
  universe_genes,
  contrast,
  gene_set,
  direction
) {
  selected_genes <- sort(unique(selected_genes))
  universe_genes <- sort(unique(universe_genes))

  selected_annotated <- intersect(selected_genes, term2gene$gene)
  universe_annotated <- intersect(universe_genes, term2gene$gene)

  universe_pairs <- term2gene[
    term2gene$gene %in% universe_annotated,
    ,
    drop = FALSE
  ]
  pathway_sizes <- table(universe_pairs$term)
  testable_pathways <- names(pathway_sizes)[
    pathway_sizes >= MIN_GENE_SET_SIZE &
      pathway_sizes <= MAX_GENE_SET_SIZE
  ]

  summary_row <- data.frame(
    contrast = contrast,
    gene_set = gene_set,
    direction = direction,
    selected_genes_total = length(selected_genes),
    selected_genes_KEGG_annotated = length(selected_annotated),
    universe_genes_total = length(universe_genes),
    universe_genes_KEGG_annotated = length(universe_annotated),
    KEGG_pathways_in_annotation = length(unique(term2gene$term)),
    KEGG_pathways_within_size_limits = length(testable_pathways),
    enriched_pathways_padj_0.05 = 0L,
    status = "not_run",
    stringsAsFactors = FALSE
  )

  if (length(selected_annotated) < MIN_GENE_SET_SIZE) {
    summary_row$status <- "too_few_annotated_selected_genes"
    return(list(result = NULL, summary = summary_row))
  }

  enricher_args <- list(
    gene = selected_annotated,
    universe = universe_annotated,
    TERM2GENE = term2gene,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = MIN_GENE_SET_SIZE,
    maxGSSize = MAX_GENE_SET_SIZE
  )

  if (!is.null(term2name)) {
    enricher_args$TERM2NAME <- term2name
  }

  x <- do.call(clusterProfiler::enricher, enricher_args)

  if (is.null(x)) {
    summary_row$status <- "no_result_object"
    return(list(result = NULL, summary = summary_row))
  }

  result <- as.data.frame(x)
  if (nrow(result) == 0L) {
    summary_row$status <- "no_pathways_tested"
    return(list(result = result, summary = summary_row))
  }

  result$contrast <- contrast
  result$gene_set <- gene_set
  result$direction <- direction
  result$significant_padj_0.05 <- (
    !is.na(result$p.adjust) & result$p.adjust < PADJ_CUTOFF
  )

  if (!is.null(metadata_unique)) {
    result <- merge(
      result,
      metadata_unique,
      by.x = "ID",
      by.y = "pathway",
      all.x = TRUE,
      sort = FALSE
    )
  }

  result <- result[
    order(result$p.adjust, result$pvalue, -result$Count),
    ,
    drop = FALSE
  ]
  rownames(result) <- NULL

  summary_row$enriched_pathways_padj_0.05 <- sum(
    result$significant_padj_0.05,
    na.rm = TRUE
  )
  summary_row$status <- "completed"

  list(result = result, summary = summary_row)
}

make_dotplot <- function(result, path, title_text, top_n = 15L) {
  if (is.null(result) || nrow(result) == 0L) {
    return(invisible(FALSE))
  }

  plot_df <- result[
    !is.na(result$p.adjust) & result$p.adjust < PADJ_CUTOFF,
    ,
    drop = FALSE
  ]
  if (nrow(plot_df) == 0L) return(invisible(FALSE))

  plot_df <- head(plot_df, top_n)
  plot_df$gene_ratio_numeric <- parse_ratio(plot_df$GeneRatio)
  plot_df$minus_log10_padj <- -log10(plot_df$p.adjust)
  plot_df$Description <- factor(
    plot_df$Description,
    levels = rev(plot_df$Description)
  )

  p <- ggplot(
    plot_df,
    aes(
      x = gene_ratio_numeric,
      y = Description,
      size = Count,
      colour = minus_log10_padj
    )
  ) +
    geom_point() +
    labs(
      title = title_text,
      x = "Gene ratio",
      y = NULL,
      size = "Gene count",
      colour = expression(-log[10]("adjusted p"))
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 8)
    )

  ggsave(
    filename = path,
    plot = p,
    width = 8.5,
    height = 6.5,
    units = "in"
  )

  invisible(TRUE)
}

# Read and validate DESeq2 contrast.
required_de_columns <- c(
  "gene_id", "log2FoldChange", "pvalue", "padj"
)
contrast_data <- lapply(
  contrast_files,
  read_tsv_checked,
  required_columns = required_de_columns
)

for (contrast in names(contrast_data)) {
  contrast_data[[contrast]]$gene_id <- trimws(
    as.character(contrast_data[[contrast]]$gene_id)
  )

  if (any(is.na(contrast_data[[contrast]]$gene_id) |
          contrast_data[[contrast]]$gene_id == "")) {
    stop_with("Missing gene IDs in DESeq2 contrast ", contrast, ".")
  }
  if (anyDuplicated(contrast_data[[contrast]]$gene_id)) {
    stop_with("Duplicated gene IDs in ", contrast, ".")
  }

  nonstructural <- setdiff(
    contrast_data[[contrast]]$gene_id,
    structural_genes
  )
  if (length(nonstructural) > 0L) {
    stop_with(
      "DESeq2 genes absent from structural annotation in ",
      contrast, ": ",
      paste(head(nonstructural, 20L), collapse = ", ")
    )
  }
}

tested_gene_sets <- lapply(
  contrast_data,
  function(x) sort(unique(x$gene_id))
)
tested_counts <- vapply(tested_gene_sets, length, integer(1))

if (any(tested_counts != EXPECTED_TESTED_GENES)) {
  stop_with(
    "Expected ", EXPECTED_TESTED_GENES,
    " tested genes per contrast, observed: ",
    paste(
      names(tested_counts),
      tested_counts,
      sep = "=",
      collapse = ", "
    )
  )
}

reference_universe <- tested_gene_sets[[1]]

if (length(tested_gene_sets) > 1L) {
  for (contrast in names(tested_gene_sets)[-1]) {
    if (!identical(reference_universe, tested_gene_sets[[contrast]])) {
      stop_with("The tested-gene universe differs between contrasts.")
    }
  }
}

selection_definitions <- list(
  significant = function(df) {
    !is.na(df$padj) & df$padj < PADJ_CUTOFF
  },
  strong = function(df) {
    !is.na(df$padj) &
      df$padj < PADJ_CUTOFF &
      !is.na(df$log2FoldChange) &
      abs(df$log2FoldChange) > STRONG_LFC_CUTOFF
  }
)

all_results <- list()
all_summaries <- list()
result_index <- 0L
summary_index <- 0L

for (contrast in names(contrast_data)) {
  df <- contrast_data[[contrast]]

  for (gene_set in names(selection_definitions)) {
    selected_flag <- selection_definitions[[gene_set]](df)

    for (direction in c("up", "down")) {
      direction_flag <- if (direction == "up") {
        !is.na(df$log2FoldChange) & df$log2FoldChange > 0
      } else {
        !is.na(df$log2FoldChange) & df$log2FoldChange < 0
      }

      selected_genes <- df$gene_id[selected_flag & direction_flag]

      run <- run_enrichment(
        selected_genes = selected_genes,
        universe_genes = reference_universe,
        contrast = contrast,
        gene_set = gene_set,
        direction = direction
      )

      summary_index <- summary_index + 1L
      all_summaries[[summary_index]] <- run$summary

      prefix <- paste(contrast, gene_set, direction, sep = "__")
      all_path <- file.path(
        output_dir, paste0(prefix, "__all_terms.tsv")
      )
      sig_path <- file.path(
        output_dir, paste0(prefix, "__significant_terms.tsv")
      )

      if (is.null(run$result)) {
        write_empty_tsv(all_path)
        write_empty_tsv(sig_path)
        next
      }

      write.table(
        run$result,
        file = all_path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )

      significant_result <- run$result[
        run$result$significant_padj_0.05,
        ,
        drop = FALSE
      ]
      write.table(
        significant_result,
        file = sig_path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )

      if (nrow(run$result) > 0L) {
        result_index <- result_index + 1L
        all_results[[result_index]] <- run$result
      }

      make_dotplot(
        run$result,
        file.path(output_dir, paste0(prefix, "__dotplot.pdf")),
        paste(
          contrast,
          paste0("(", gene_set, ", ", direction, ")"),
          "KEGG enrichment"
        )
      )
    }
  }
}

summary_df <- do.call(rbind, all_summaries)
write.table(
  summary_df,
  file = file.path(
    output_dir, "kegg_enrichment_run_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

combined_results <- if (length(all_results) > 0L) {
  do.call(rbind, all_results)
} else {
  data.frame()
}
rownames(combined_results) <- NULL

write.table(
  combined_results,
  file = file.path(
    output_dir, "all_KEGG_enrichment_results.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

significant_combined <- if (nrow(combined_results) > 0L) {
  combined_results[
    !is.na(combined_results$p.adjust) &
      combined_results$p.adjust < PADJ_CUTOFF,
    ,
    drop = FALSE
  ]
} else {
  data.frame()
}
rownames(significant_combined) <- NULL

write.table(
  significant_combined,
  file = file.path(
    output_dir, "all_significant_KEGG_pathways.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

universe_annotated <- intersect(
  reference_universe,
  term2gene$gene
)

identical_universe <- if (length(tested_gene_sets) > 1L) {
  all(vapply(
    tested_gene_sets,
    identical,
    logical(1),
    reference_universe
  ))
} else {
  TRUE
}

validation <- data.frame(
  check = c(
    "structural_genes_in_annotation_table",
    "contrasts_present",
    "tested_genes_per_contrast",
    "identical_universe_across_contrasts",
    "all_tested_genes_present_in_structural_annotation",
    "KEGG_annotated_genes_in_full_TERM2GENE",
    "gene_pathway_pairs_in_full_TERM2GENE",
    "pathways_in_full_TERM2GENE",
    "KEGG_annotated_genes_in_tested_universe",
    "pathway_metadata_available",
    "TERM2NAME_used",
    "primary_padj_cutoff",
    "strong_abs_log2fc_cutoff",
    "min_gene_set_size",
    "max_gene_set_size"
  ),
  value = c(
    nrow(gene_annotation),
    length(contrast_data),
    paste(tested_counts, collapse = ";"),
    identical_universe,
    all(reference_universe %in% structural_genes),
    length(unique(term2gene$gene)),
    nrow(term2gene),
    length(unique(term2gene$term)),
    length(universe_annotated),
    metadata_available,
    term2name_used,
    PADJ_CUTOFF,
    STRONG_LFC_CUTOFF,
    MIN_GENE_SET_SIZE,
    MAX_GENE_SET_SIZE
  ),
  stringsAsFactors = FALSE
)

write.table(
  validation,
  file = file.path(
    output_dir,
    "22_run_fusarium_kegg_enrichment.validation.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    output_dir,
    "22_run_fusarium_kegg_enrichment.sessionInfo.txt"
  )
)

message("")
message("Fusarium KEGG enrichment completed successfully.")
message("Tested-gene universe:       ", length(reference_universe))
message("KEGG-annotated universe:    ", length(universe_annotated))
message("Full KEGG pathway genes:    ", length(unique(term2gene$gene)))
message("Full gene-pathway pairs:    ", nrow(term2gene))
message("Full KEGG pathways:         ", length(unique(term2gene$term)))
message("Pathway metadata available: ", metadata_available)
message("TERM2NAME used:             ", term2name_used)
message("Significant KEGG rows:      ", nrow(significant_combined))
message("Results directory:          ", output_dir)
message("")
message("Primary interpretation files:")
message("  all_significant_KEGG_pathways.tsv")
message("  kegg_enrichment_run_summary.tsv")
