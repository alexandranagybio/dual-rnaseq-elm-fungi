#!/usr/bin/env Rscript

# 20_run_fusarium_go_enrichment.R
#
# GO over-representation analysis for Fusarium salinense Trinity genes.
# Original protein-level eggNOG GO annotations are mapped to validated Trinity
# genes and deduplicated before enrichment.
#
# Primary: raw_padj < 0.05, split up/down.
# Sensitivity: raw_padj < 0.05 and abs(raw_log2FoldChange) > 1.
# Universe: all 10,414 genes present in the DESeq2 dataset.

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ggplot2)
  library(AnnotationDbi)
  library(GO.db)
})

EXPECTED_STRUCTURAL_GENES <- 15192L
EXPECTED_PROTEINS <- 35327L
EXPECTED_DESEQ2_GENES <- 10414L
EXPECTED_GENES_WITH_PADJ <- 9000L
EXPECTED_GENES_WITH_NA_PADJ <- 1414L
EXPECTED_SIGNIFICANT_GENES <- 2973L
EXPECTED_SIGNIFICANT_UP <- 1550L
EXPECTED_SIGNIFICANT_DOWN <- 1423L
EXPECTED_STRONG_GENES <- 149L

PADJ_CUTOFF <- 0.05
STRONG_LFC_CUTOFF <- 1
MIN_GENE_SET_SIZE <- 5L
MAX_GENE_SET_SIZE <- 500L

functional_dir <- file.path("results", "fusarium", "functional_annotation")
publication_dir <- file.path("results", "fusarium", "publication_annotation")
output_dir <- file.path("results", "fusarium", "go_enrichment")
tables_dir <- file.path(output_dir, "tables")
plots_dir <- file.path(output_dir, "plots")
diagnostics_dir <- file.path(output_dir, "diagnostics")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)

protein_file <- file.path(
  functional_dir, "tables", "fusarium_protein_functional_annotation.tsv"
)
gene_file <- file.path(
  functional_dir, "tables", "fusarium_gene_functional_annotation.tsv"
)
publication_file <- file.path(
  publication_dir, "tables", "fusarium_publication_annotation_full.tsv"
)
eggnog_file <- file.path(
  "data", "external", "fusarium_annotation",
  "fusarium_complete.emapper.annotations"
)

stop_with <- function(...) stop(paste0(...), call. = FALSE)

read_tsv_checked <- function(path, required = NULL) {
  if (!file.exists(path)) stop_with("Missing input file: ", path)
  x <- read.delim(
    path, sep = "\t", header = TRUE, quote = "", comment.char = "",
    check.names = FALSE, na.strings = c("", "NA", "NaN")
  )
  missing <- setdiff(required, colnames(x))
  if (length(missing)) {
    stop_with("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  }
  x
}

as_logical_checked <- function(x, field) {
  if (is.logical(x)) return(x)
  y <- toupper(trimws(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("TRUE", "T", "1")] <- TRUE
  out[y %in% c("FALSE", "F", "0")] <- FALSE
  bad <- is.na(out) & !is.na(x) & trimws(as.character(x)) != ""
  if (any(bad)) stop_with("Invalid logical values in ", field)
  out
}

read_eggnog <- function(path) {
  if (!file.exists(path)) stop_with("Missing eggNOG file: ", path)
  header_lines <- grep("^#query\t", readLines(path, warn = FALSE), value = TRUE)
  if (length(header_lines) != 1L) {
    stop_with("Expected one #query header; observed ", length(header_lines))
  }
  header <- strsplit(sub("^#", "", header_lines), "\t")[[1]]
  x <- read.delim(
    path, sep = "\t", header = FALSE, comment.char = "#", quote = "",
    fill = TRUE, check.names = FALSE, na.strings = c("", "-", "NA")
  )
  if (ncol(x) != length(header)) stop_with("eggNOG column-count mismatch")
  colnames(x) <- header
  if (!all(c("query", "GOs") %in% colnames(x))) {
    stop_with("eggNOG file lacks #query and/or GOs")
  }
  x
}

split_go_rows <- function(eggnog) {
  x <- eggnog[!is.na(eggnog$GOs) & trimws(eggnog$GOs) != "", c("query", "GOs")]
  if (!nrow(x)) stop_with("No GO annotations found in eggNOG file")
  pieces <- strsplit(x$GOs, ",", fixed = TRUE)
  out <- data.frame(
    protein_id = rep(x[["query"]], lengths(pieces)),
    go_id = trimws(unlist(pieces, use.names = FALSE)),
    stringsAsFactors = FALSE
  )
  unique(out[grepl("^GO:[0-9]{7}$", out$go_id), ])
}

resolve_go <- function(go_ids) {
  go_ids <- sort(unique(go_ids))
  metadata <- AnnotationDbi::select(
    GO.db, keys = go_ids, keytype = "GOID", columns = c("TERM", "ONTOLOGY")
  )
  colnames(metadata) <- c("go_id", "go_name", "ontology")
  metadata <- unique(metadata)
  resolved <- metadata[
    !is.na(metadata$ontology) & metadata$ontology %in% c("BP", "MF", "CC"),
  ]
  unresolved <- data.frame(
    go_id = setdiff(go_ids, resolved$go_id), stringsAsFactors = FALSE
  )
  list(resolved = resolved, unresolved = unresolved)
}

protein_map <- read_tsv_checked(protein_file, c("protein_id", "gene_id"))
gene_map <- read_tsv_checked(gene_file, "gene_id")

if (nrow(protein_map) != EXPECTED_PROTEINS) {
  stop_with("Expected ", EXPECTED_PROTEINS, " proteins; observed ", nrow(protein_map))
}
if (anyDuplicated(protein_map$protein_id)) stop_with("Duplicated protein IDs")
if (nrow(gene_map) != EXPECTED_STRUCTURAL_GENES) {
  stop_with("Expected ", EXPECTED_STRUCTURAL_GENES, " genes; observed ", nrow(gene_map))
}
if (anyDuplicated(gene_map$gene_id)) stop_with("Duplicated gene IDs")

eggnog <- read_eggnog(eggnog_file)
if (anyDuplicated(eggnog[["#query"]])) stop_with("Duplicated eggNOG query IDs")

protein_go <- split_go_rows(eggnog)
missing_proteins <- setdiff(unique(protein_go$protein_id), protein_map$protein_id)
if (length(missing_proteins)) {
  stop_with(length(missing_proteins), " GO-annotated proteins absent from proteome")
}

protein_go <- merge(
  protein_go, protein_map[, c("protein_id", "gene_id")],
  by = "protein_id", all.x = TRUE, sort = FALSE
)
if (anyNA(protein_go$gene_id)) stop_with("Unmapped GO-annotated proteins")

go_metadata <- resolve_go(protein_go$go_id)
write.table(
  go_metadata$unresolved,
  file.path(diagnostics_dir, "unresolved_or_obsolete_GO_ids.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

go_long <- merge(
  protein_go, go_metadata$resolved,
  by = "go_id", all = FALSE, sort = FALSE
)
go_long <- unique(go_long[, c("gene_id", "protein_id", "go_id", "go_name", "ontology")])
term2gene_long <- unique(go_long[, c("gene_id", "go_id", "go_name", "ontology")])

write.table(
  go_long, file.path(tables_dir, "fusarium_GO_protein_to_gene_long.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  term2gene_long, file.path(tables_dir, "fusarium_GO_TERM2GENE_long.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

de <- read_tsv_checked(
  publication_file,
  c("gene_id", "in_deseq2_dataset", "raw_log2FoldChange", "raw_padj")
)
if (nrow(de) != EXPECTED_STRUCTURAL_GENES) stop_with("Unexpected publication rows")
if (anyDuplicated(de$gene_id)) stop_with("Duplicated publication gene IDs")

de$in_deseq2_dataset <- as_logical_checked(de$in_deseq2_dataset, "in_deseq2_dataset")
tested <- de[!is.na(de$in_deseq2_dataset) & de$in_deseq2_dataset, ]
if (nrow(tested) != EXPECTED_DESEQ2_GENES) {
  stop_with("Expected ", EXPECTED_DESEQ2_GENES, " tested genes; observed ", nrow(tested))
}

significant <- !is.na(tested$raw_padj) & tested$raw_padj < PADJ_CUTOFF
strong <- significant & !is.na(tested$raw_log2FoldChange) &
  abs(tested$raw_log2FoldChange) > STRONG_LFC_CUTOFF

observed <- c(
  genes_with_padj = sum(!is.na(tested$raw_padj)),
  genes_with_na_padj = sum(is.na(tested$raw_padj)),
  significant = sum(significant),
  significant_up = sum(significant & tested$raw_log2FoldChange > 0, na.rm = TRUE),
  significant_down = sum(significant & tested$raw_log2FoldChange < 0, na.rm = TRUE),
  strong = sum(strong)
)
expected <- c(
  genes_with_padj = EXPECTED_GENES_WITH_PADJ,
  genes_with_na_padj = EXPECTED_GENES_WITH_NA_PADJ,
  significant = EXPECTED_SIGNIFICANT_GENES,
  significant_up = EXPECTED_SIGNIFICANT_UP,
  significant_down = EXPECTED_SIGNIFICANT_DOWN,
  strong = EXPECTED_STRONG_GENES
)
if (any(observed != expected)) {
  bad <- names(observed)[observed != expected]
  stop_with(
    "DESeq2 validation failed: ",
    paste(paste0(bad, "=", observed[bad], " expected ", expected[bad]), collapse = "; ")
  )
}

run_enrichment <- function(selected_genes, universe_genes, ontology, set_name, direction) {
  m <- term2gene_long[term2gene_long$ontology == ontology, ]
  term2gene <- unique(m[, c("go_id", "gene_id")])
  term2name <- unique(m[, c("go_id", "go_name")])
  colnames(term2gene) <- c("term", "gene")
  colnames(term2name) <- c("term", "name")

  selected_annotated <- intersect(unique(selected_genes), term2gene$gene)
  universe_annotated <- intersect(unique(universe_genes), term2gene$gene)

  metadata <- data.frame(
    contrast = "interaction_vs_self", gene_set = set_name,
    direction = direction, ontology = ontology,
    selected_genes_total = length(unique(selected_genes)),
    selected_genes_GO_annotated = length(selected_annotated),
    universe_genes_total = length(unique(universe_genes)),
    universe_genes_GO_annotated = length(universe_annotated),
    available_GO_terms = length(unique(term2gene$term)),
    tested_GO_terms = 0L, enriched_terms_padj_0.05 = 0L,
    status = "not_run", stringsAsFactors = FALSE
  )

  if (length(selected_annotated) < MIN_GENE_SET_SIZE) {
    metadata$status <- "too_few_annotated_selected_genes"
    return(list(result = NULL, metadata = metadata))
  }

  enrichment <- enricher(
    gene = selected_annotated, universe = universe_annotated,
    TERM2GENE = term2gene, TERM2NAME = term2name,
    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
    minGSSize = MIN_GENE_SET_SIZE, maxGSSize = MAX_GENE_SET_SIZE
  )

  if (is.null(enrichment)) {
    metadata$status <- "no_result_object"
    return(list(result = NULL, metadata = metadata))
  }

  result <- as.data.frame(enrichment)
  if (!nrow(result)) {
    metadata$status <- "no_terms_tested"
    return(list(result = result, metadata = metadata))
  }

  result$contrast <- "interaction_vs_self"
  result$gene_set <- set_name
  result$direction <- direction
  result$ontology <- ontology
  result$significant_padj_0.05 <- !is.na(result$p.adjust) & result$p.adjust < PADJ_CUTOFF
  result <- result[order(result$p.adjust, result$pvalue, -result$Count), ]

  metadata$tested_GO_terms <- nrow(result)
  metadata$enriched_terms_padj_0.05 <- sum(result$significant_padj_0.05)
  metadata$status <- "completed"
  list(result = result, metadata = metadata)
}

ratio_numeric <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(z) {
    if (length(z) != 2L) return(NA_real_)
    as.numeric(z[1]) / as.numeric(z[2])
  }, numeric(1))
}

make_dotplot <- function(result, path, title, top_n = 15L) {
  if (is.null(result) || !nrow(result)) return(invisible(FALSE))
  d <- result[!is.na(result$p.adjust) & result$p.adjust < PADJ_CUTOFF, ]
  if (!nrow(d)) return(invisible(FALSE))
  d <- head(d, top_n)
  d$gene_ratio_numeric <- ratio_numeric(d$GeneRatio)
  d$minus_log10_padj <- -log10(d$p.adjust)
  d$Description <- factor(d$Description, levels = rev(d$Description))
  p <- ggplot(d, aes(gene_ratio_numeric, Description, size = Count, colour = minus_log10_padj)) +
    geom_point() +
    labs(
      title = title, x = "Gene ratio", y = NULL,
      size = "Gene count", colour = expression(-log[10]("adjusted p"))
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"), axis.text.y = element_text(size = 8))
  ggsave(path, p, width = 8.5, height = 6.5, units = "in")
  invisible(TRUE)
}

universe <- sort(unique(tested$gene_id))
sets <- list(significant = significant, strong = strong)
all_results <- list()
all_metadata <- list()
ri <- 0L
mi <- 0L

for (set_name in names(sets)) {
  for (direction in c("up", "down")) {
    direction_flag <- if (direction == "up") {
      !is.na(tested$raw_log2FoldChange) & tested$raw_log2FoldChange > 0
    } else {
      !is.na(tested$raw_log2FoldChange) & tested$raw_log2FoldChange < 0
    }
    selected_genes <- tested$gene_id[sets[[set_name]] & direction_flag]

    for (ontology in c("BP", "MF", "CC")) {
      run <- run_enrichment(selected_genes, universe, ontology, set_name, direction)
      mi <- mi + 1L
      all_metadata[[mi]] <- run$metadata
      prefix <- paste("interaction_vs_self", set_name, direction, ontology, sep = "__")
      all_path <- file.path(tables_dir, paste0(prefix, "__all_terms.tsv"))
      sig_path <- file.path(tables_dir, paste0(prefix, "__significant_terms.tsv"))

      if (is.null(run$result)) {
        write.table(data.frame(), all_path, sep = "\t", quote = FALSE, row.names = FALSE)
        write.table(data.frame(), sig_path, sep = "\t", quote = FALSE, row.names = FALSE)
        next
      }

      write.table(run$result, all_path, sep = "\t", quote = FALSE, row.names = FALSE)
      sig_result <- run$result[run$result$significant_padj_0.05, ]
      write.table(sig_result, sig_path, sep = "\t", quote = FALSE, row.names = FALSE)
      if (nrow(run$result)) {
        ri <- ri + 1L
        all_results[[ri]] <- run$result
      }
      make_dotplot(
        run$result, file.path(plots_dir, paste0(prefix, "__dotplot.pdf")),
        paste("Fusarium interaction vs self", paste0("(", set_name, ", ", direction, ")"), ontology, "GO enrichment")
      )
    }
  }
}

metadata_df <- do.call(rbind, all_metadata)
write.table(
  metadata_df, file.path(diagnostics_dir, "go_enrichment_run_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

combined <- if (length(all_results)) do.call(rbind, all_results) else data.frame()
write.table(
  combined, file.path(tables_dir, "all_GO_enrichment_results.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
sig_combined <- if (nrow(combined)) {
  combined[!is.na(combined$p.adjust) & combined$p.adjust < PADJ_CUTOFF, ]
} else data.frame()
write.table(
  sig_combined, file.path(tables_dir, "all_significant_GO_terms.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

validation <- data.frame(
  check = c(
    "structural_gene_rows", "protein_rows", "deseq2_gene_rows",
    "genes_with_padj", "genes_with_na_padj", "significant_padj_lt_0.05",
    "significant_up", "significant_down", "strong_raw_abs_lfc_gt_1",
    "eggnog_annotation_rows", "proteins_with_GO",
    "unique_GO_ids_before_GOdb_resolution", "unresolved_or_obsolete_GO_ids",
    "GO_annotated_proteins_after_resolution", "GO_annotated_structural_genes",
    "GO_annotated_DESeq2_universe_genes"
  ),
  value = c(
    nrow(gene_map), nrow(protein_map), nrow(tested), observed,
    nrow(eggnog), length(unique(protein_go$protein_id)),
    length(unique(protein_go$go_id)), nrow(go_metadata$unresolved),
    length(unique(go_long$protein_id)), length(unique(term2gene_long$gene_id)),
    length(intersect(universe, term2gene_long$gene_id))
  ),
  expected = c(
    EXPECTED_STRUCTURAL_GENES, EXPECTED_PROTEINS, EXPECTED_DESEQ2_GENES,
    expected, rep(NA, 7)
  ),
  stringsAsFactors = FALSE
)
validation$status <- ifelse(
  is.na(validation$expected), "INFO",
  ifelse(as.character(validation$value) == as.character(validation$expected), "PASS", "FAIL")
)
write.table(
  validation,
  file.path(diagnostics_dir, "20_run_fusarium_go_enrichment.validation.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

run_info <- data.frame(
  field = c(
    "script", "run_timestamp_utc", "contrast", "analysis_level", "GO_source",
    "GO_metadata_source", "background_definition", "primary_selection",
    "sensitivity_selection", "direction_definition", "ontology_handling",
    "multiple_testing", "protein_annotation_file", "gene_annotation_file",
    "publication_file", "eggnog_file"
  ),
  value = c(
    "20_run_fusarium_go_enrichment.R",
    format(Sys.time(), tz = "UTC", usetz = TRUE),
    "interaction_vs_self", "Trinity gene",
    "original protein-level eggNOG-mapper annotations", "GO.db",
    "all 10414 genes present in DESeq2, restricted per ontology to GO-annotated genes",
    "raw_padj < 0.05",
    "raw_padj < 0.05 and abs(raw_log2FoldChange) > 1",
    "sign of raw_log2FoldChange", "BP, MF, and CC separately",
    "Benjamini-Hochberg", normalizePath(protein_file), normalizePath(gene_file),
    normalizePath(publication_file), normalizePath(eggnog_file)
  ),
  stringsAsFactors = FALSE
)
write.table(run_info, file.path(output_dir, "run_info.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(
  capture.output(sessionInfo()),
  file.path(diagnostics_dir, "20_run_fusarium_go_enrichment.sessionInfo.txt")
)

if (any(validation$status == "FAIL")) {
  stop_with("Validation failed; inspect diagnostics/20_run_fusarium_go_enrichment.validation.tsv")
}

message("")
message("Fusarium GO enrichment completed successfully.")
message("DESeq2 universe:       ", length(universe))
message("GO-annotated universe: ", length(intersect(universe, term2gene_long$gene_id)))
message("GO-annotated genes:    ", length(unique(term2gene_long$gene_id)))
message("Unresolved GO IDs:     ", nrow(go_metadata$unresolved))
message("Significant GO rows:   ", nrow(sig_combined))
message("Validation status:     PASS")
message("Results directory:     ", output_dir)
message("")
message("Primary interpretation files:")
message("  tables/all_significant_GO_terms.tsv")
message("  diagnostics/go_enrichment_run_summary.tsv")
