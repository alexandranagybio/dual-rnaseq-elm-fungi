#!/usr/bin/env Rscript

# 10_ophiostoma_go_enrichment.R
#
# Gene Ontology over-representation analysis for rebuilt Ophiostoma novo-ulmi
# gene-level DESeq2 results.
#
# This corrected version builds TERM2GENE and TERM2NAME directly from the
# original JGI GO table rather than from collapsed annotation columns.
#
# Primary analysis:
#   - significant genes: padj < 0.05
#   - split by log2FoldChange direction
#   - background: all 8,560 genes tested by DESeq2
#   - GO enrichment performed separately for BP, MF, and CC
#
# Sensitivity analysis:
#   - strong genes: padj < 0.05 and abs(log2FoldChange) > 1

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ggplot2)
})

EXPECTED_TESTED_GENES <- 8560L
EXPECTED_STRUCTURAL_GENES <- 8640L
EXPECTED_GO_ANNOTATED_GENES <- 4680L

PADJ_CUTOFF <- 0.05
STRONG_LFC_CUTOFF <- 1
MIN_GENE_SET_SIZE <- 5L
MAX_GENE_SET_SIZE <- 500L

annotation_dir <- file.path("results", "ophiostoma", "annotation")
output_dir <- file.path("results", "ophiostoma", "go_enrichment")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gene_map_file <- file.path(
  annotation_dir,
  "ophiostoma_gene_annotation_map.tsv"
)

go_file <- file.path(
  "data",
  "annotation",
  "ophiostoma",
  "GO",
  "Ophnu1_GeneCatalog_proteins_20170425_GO.tab"
)

contrast_files <- c(
  interaction_vs_self = file.path(
    annotation_dir,
    "interaction_vs_self_all_genes.annotated.tsv"
  ),
  onu_vs_self = file.path(
    annotation_dir,
    "onu_vs_self_all_genes.annotated.tsv"
  ),
  interaction_vs_onu = file.path(
    annotation_dir,
    "interaction_vs_onu_all_genes.annotated.tsv"
  )
)

stop_with <- function(...) {
  stop(paste0(...), call. = FALSE)
}

read_tsv_checked <- function(path, required_columns = NULL) {
  if (!file.exists(path)) {
    stop_with("Missing input file: ", path)
  }

  x <- read.delim(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    na.strings = c("", "NA", "NaN")
  )

  if (!is.null(required_columns)) {
    missing_columns <- setdiff(required_columns, colnames(x))
    if (length(missing_columns) > 0L) {
      stop_with(
        "Input ", path, " is missing required columns: ",
        paste(missing_columns, collapse = ", ")
      )
    }
  }

  x
}

normalize_ontology <- function(x) {
  y <- tolower(trimws(x))
  out <- rep(NA_character_, length(y))
  out[y %in% c("biological_process", "biological process", "bp")] <- "BP"
  out[y %in% c("molecular_function", "molecular function", "mf")] <- "MF"
  out[y %in% c("cellular_component", "cellular component", "cc")] <- "CC"
  out
}

build_go_mapping <- function(gene_map_file, go_file) {
  gene_map <- read_tsv_checked(
    gene_map_file,
    required_columns = c("gene_id", "protein_id")
  )

  if (nrow(gene_map) != EXPECTED_STRUCTURAL_GENES) {
    stop_with(
      "Expected ", EXPECTED_STRUCTURAL_GENES,
      " structural genes in annotation map, observed ", nrow(gene_map), "."
    )
  }

  if (anyDuplicated(gene_map$gene_id)) {
    stop_with("Duplicated gene_id values in annotation map.")
  }

  if (anyDuplicated(gene_map$protein_id)) {
    stop_with("Duplicated protein_id values in annotation map.")
  }

  go <- read.delim(
    go_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  colnames(go)[1] <- sub("^#", "", colnames(go)[1])

  required_go_columns <- c(
    "proteinId",
    "gotermId",
    "goName",
    "gotermType",
    "goAcc"
  )

  missing_go_columns <- setdiff(required_go_columns, colnames(go))
  if (length(missing_go_columns) > 0L) {
    stop_with(
      "GO table is missing required columns: ",
      paste(missing_go_columns, collapse = ", ")
    )
  }

  go$proteinId <- as.character(go$proteinId)
  gene_map$protein_id <- as.character(gene_map$protein_id)

  go_long <- merge(
    go,
    gene_map[, c("gene_id", "protein_id")],
    by.x = "proteinId",
    by.y = "protein_id",
    all.x = FALSE,
    all.y = FALSE,
    sort = FALSE
  )

  go_long$ontology <- normalize_ontology(go_long$gotermType)

  invalid_ontology <- is.na(go_long$ontology)
  if (any(invalid_ontology)) {
    bad_types <- unique(go_long$gotermType[invalid_ontology])
    stop_with(
      "Unrecognized GO ontology labels: ",
      paste(bad_types, collapse = ", ")
    )
  }

  go_long <- unique(
    go_long[, c(
      "gene_id",
      "proteinId",
      "goAcc",
      "goName",
      "gotermType",
      "ontology"
    )]
  )

  colnames(go_long) <- c(
    "gene_id",
    "protein_id",
    "go_id",
    "go_name",
    "go_type",
    "ontology"
  )

  annotated_gene_count <- length(unique(go_long$gene_id))
  if (annotated_gene_count != EXPECTED_GO_ANNOTATED_GENES) {
    warning(
      "Expected ", EXPECTED_GO_ANNOTATED_GENES,
      " GO-annotated genes, observed ", annotated_gene_count, "."
    )
  }

  list(
    gene_map = gene_map,
    go_long = go_long
  )
}

run_enrichment <- function(
  selected_genes,
  universe_genes,
  go_long,
  ontology,
  contrast,
  set_name,
  direction
) {
  ontology_map <- go_long[go_long$ontology == ontology, , drop = FALSE]

  term2gene <- unique(ontology_map[, c("go_id", "gene_id")])
  term2name <- unique(ontology_map[, c("go_id", "go_name")])

  colnames(term2gene) <- c("term", "gene")
  colnames(term2name) <- c("term", "name")

  selected_annotated <- intersect(selected_genes, term2gene$gene)
  universe_annotated <- intersect(universe_genes, term2gene$gene)

  metadata <- data.frame(
    contrast = contrast,
    gene_set = set_name,
    direction = direction,
    ontology = ontology,
    selected_genes_total = length(selected_genes),
    selected_genes_GO_annotated = length(selected_annotated),
    universe_genes_total = length(universe_genes),
    universe_genes_GO_annotated = length(universe_annotated),
    tested_GO_terms = length(unique(term2gene$term)),
    enriched_terms_padj_0.05 = 0L,
    status = "not_run",
    stringsAsFactors = FALSE
  )

  if (length(selected_annotated) < MIN_GENE_SET_SIZE) {
    metadata$status <- "too_few_annotated_selected_genes"
    return(list(result = NULL, metadata = metadata))
  }

  enrichment <- enricher(
    gene = selected_annotated,
    universe = universe_annotated,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = MIN_GENE_SET_SIZE,
    maxGSSize = MAX_GENE_SET_SIZE
  )

  if (is.null(enrichment)) {
    metadata$status <- "no_result_object"
    return(list(result = NULL, metadata = metadata))
  }

  result <- as.data.frame(enrichment)

  if (nrow(result) == 0L) {
    metadata$status <- "no_terms_tested"
    return(list(result = result, metadata = metadata))
  }

  result$contrast <- contrast
  result$gene_set <- set_name
  result$direction <- direction
  result$ontology <- ontology
  result$significant_padj_0.05 <- !is.na(result$p.adjust) &
    result$p.adjust < PADJ_CUTOFF

  result <- result[
    order(result$p.adjust, result$pvalue, -result$Count),
    ,
    drop = FALSE
  ]

  metadata$enriched_terms_padj_0.05 <- sum(
    result$significant_padj_0.05,
    na.rm = TRUE
  )
  metadata$status <- "completed"

  list(result = result, metadata = metadata)
}

make_dotplot <- function(result_df, output_pdf, title_text, top_n = 15L) {
  if (is.null(result_df) || nrow(result_df) == 0L) {
    return(invisible(FALSE))
  }

  significant <- result_df[
    !is.na(result_df$p.adjust) & result_df$p.adjust < PADJ_CUTOFF,
    ,
    drop = FALSE
  ]

  if (nrow(significant) == 0L) {
    return(invisible(FALSE))
  }

  plot_df <- head(significant, top_n)
  plot_df$minus_log10_padj <- -log10(plot_df$p.adjust)
  plot_df$Description <- factor(
    plot_df$Description,
    levels = rev(plot_df$Description)
  )

  p <- ggplot(
    plot_df,
    aes(
      x = GeneRatio,
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
    filename = output_pdf,
    plot = p,
    width = 8.5,
    height = 6.5,
    units = "in"
  )

  invisible(TRUE)
}

mapping <- build_go_mapping(gene_map_file, go_file)
go_long <- mapping$go_long

write.table(
  go_long,
  file = file.path(output_dir, "ophiostoma_GO_TERM2GENE_long.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

required_de_columns <- c(
  "gene_id",
  "log2FoldChange",
  "pvalue",
  "padj"
)

all_contrast_data <- lapply(
  contrast_files,
  read_tsv_checked,
  required_columns = required_de_columns
)

for (contrast in names(all_contrast_data)) {
  df <- all_contrast_data[[contrast]]
  if (anyDuplicated(df$gene_id)) {
    stop_with("Duplicated gene IDs in ", contrast, ".")
  }
}

tested_gene_sets <- lapply(
  all_contrast_data,
  function(x) sort(unique(x$gene_id))
)

tested_counts <- vapply(tested_gene_sets, length, integer(1))

if (any(tested_counts != EXPECTED_TESTED_GENES)) {
  stop_with(
    "Expected ", EXPECTED_TESTED_GENES,
    " tested genes per contrast, observed: ",
    paste(names(tested_counts), tested_counts, sep = "=", collapse = ", ")
  )
}

reference_universe <- tested_gene_sets[[1]]

for (contrast in names(tested_gene_sets)[-1]) {
  if (!identical(reference_universe, tested_gene_sets[[contrast]])) {
    stop_with("The tested-gene universe differs between contrasts.")
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
all_metadata <- list()
result_index <- 0L
metadata_index <- 0L

for (contrast in names(all_contrast_data)) {
  df <- all_contrast_data[[contrast]]

  for (set_name in names(selection_definitions)) {
    selected_flag <- selection_definitions[[set_name]](df)

    for (direction in c("up", "down")) {
      direction_flag <- if (direction == "up") {
        !is.na(df$log2FoldChange) & df$log2FoldChange > 0
      } else {
        !is.na(df$log2FoldChange) & df$log2FoldChange < 0
      }

      selected_genes <- df$gene_id[selected_flag & direction_flag]

      for (ontology in c("BP", "MF", "CC")) {
        run <- run_enrichment(
          selected_genes = selected_genes,
          universe_genes = reference_universe,
          go_long = go_long,
          ontology = ontology,
          contrast = contrast,
          set_name = set_name,
          direction = direction
        )

        metadata_index <- metadata_index + 1L
        all_metadata[[metadata_index]] <- run$metadata

        file_prefix <- paste(
          contrast,
          set_name,
          direction,
          ontology,
          sep = "__"
        )

        result_path <- file.path(
          output_dir,
          paste0(file_prefix, "__all_terms.tsv")
        )

        significant_path <- file.path(
          output_dir,
          paste0(file_prefix, "__significant_terms.tsv")
        )

        if (is.null(run$result)) {
          write.table(
            data.frame(),
            file = result_path,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
          )
          write.table(
            data.frame(),
            file = significant_path,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
          )
          next
        }

        write.table(
          run$result,
          file = result_path,
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
          file = significant_path,
          sep = "\t",
          quote = FALSE,
          row.names = FALSE
        )

        if (nrow(run$result) > 0L) {
          result_index <- result_index + 1L
          all_results[[result_index]] <- run$result
        }

        make_dotplot(
          result_df = run$result,
          output_pdf = file.path(
            output_dir,
            paste0(file_prefix, "__dotplot.pdf")
          ),
          title_text = paste(
            contrast,
            paste0("(", set_name, ", ", direction, ")"),
            ontology,
            "GO enrichment"
          )
        )
      }
    }
  }
}

metadata_df <- do.call(rbind, all_metadata)

write.table(
  metadata_df,
  file = file.path(output_dir, "go_enrichment_run_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

combined_results <- if (length(all_results) > 0L) {
  do.call(rbind, all_results)
} else {
  data.frame()
}

write.table(
  combined_results,
  file = file.path(output_dir, "all_GO_enrichment_results.tsv"),
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

write.table(
  significant_combined,
  file = file.path(output_dir, "all_significant_GO_terms.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

validation <- data.frame(
  check = c(
    "structural_genes_in_annotation_map",
    "contrasts_present",
    "tested_genes_per_contrast",
    "identical_universe_across_contrasts",
    "GO_annotated_genes_in_TERM2GENE",
    "primary_padj_cutoff",
    "strong_abs_log2fc_cutoff",
    "min_gene_set_size",
    "max_gene_set_size"
  ),
  value = c(
    nrow(mapping$gene_map),
    length(all_contrast_data),
    paste(tested_counts, collapse = ";"),
    TRUE,
    length(unique(go_long$gene_id)),
    PADJ_CUTOFF,
    STRONG_LFC_CUTOFF,
    MIN_GENE_SET_SIZE,
    MAX_GENE_SET_SIZE
  ),
  stringsAsFactors = FALSE
)

write.table(
  validation,
  file = file.path(output_dir, "10_ophiostoma_go_enrichment.validation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "10_ophiostoma_go_enrichment.sessionInfo.txt")
)

message("")
message("GO enrichment completed successfully.")
message("Tested-gene universe: ", length(reference_universe))
message("GO-annotated genes:   ", length(unique(go_long$gene_id)))
message("Significant GO rows:   ", nrow(significant_combined))
message("Results directory:     ", output_dir)
message("")
message("Primary interpretation files:")
message("  all_significant_GO_terms.tsv")
message("  go_enrichment_run_summary.tsv")
