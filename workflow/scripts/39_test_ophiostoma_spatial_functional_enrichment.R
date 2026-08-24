#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

# ==========================================================================
# Configuration
# ==========================================================================

alpha <- 0.05

classification_file <- file.path(
  "results",
  "publication",
  "ophiostoma_spatial_heatmap",
  "ophiostoma_spatial_response_classification_all.tsv"
)

annotation_file <- file.path(
  "results",
  "ophiostoma",
  "publication_annotation",
  "tables",
  "ophiostoma_publication_annotation_wide.tsv"
)

output_dir <- file.path(
  "results",
  "publication",
  "ophiostoma_spatial_functional_enrichment"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

summary_out <- file.path(
  output_dir,
  "spatial_program_functional_summary.tsv"
)

binary_out <- file.path(
  output_dir,
  "spatial_secretome_cazyme_enrichment.tsv"
)

cog_complete_out <- file.path(
  output_dir,
  "spatial_cog_enrichment_complete.tsv"
)

cog_significant_out <- file.path(
  output_dir,
  "spatial_cog_enrichment_significant.tsv"
)

run_info_out <- file.path(
  output_dir,
  "run_info.tsv"
)

# ==========================================================================
# Validation
# ==========================================================================

required_files <- c(
  classification_file,
  annotation_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required input file(s):\n",
    paste0("  - ", missing_files, collapse = "\n")
  )
}

# ==========================================================================
# Helpers
# ==========================================================================

truthy <- function(x) {
  tolower(as.character(x)) %in%
    c("true", "t", "1", "yes", "y")
}

valid_codes <- c(
  "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L",
  "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "Y", "Z"
)

expected_groups <- c(
  "Reaction-zone specific",
  "Plate-wide confrontation response",
  "Complex spatial response",
  "Non-contact-region specific"
)

# ==========================================================================
# Read spatial classification
# ==========================================================================

spatial <- read_tsv(
  classification_file,
  show_col_types = FALSE,
  progress = FALSE
) |>
  filter(!is.na(response_group)) |>
  mutate(
    response_group = as.character(response_group)
  )

if (anyDuplicated(spatial$gene_id)) {
  stop("Spatial classification contains duplicated gene IDs.")
}

observed_groups <- sort(unique(spatial$response_group))
expected_groups_sorted <- sort(expected_groups)

if (!identical(observed_groups, expected_groups_sorted)) {
  stop(
    "Unexpected spatial response groups.\nObserved: ",
    paste(observed_groups, collapse = ", ")
  )
}

# ==========================================================================
# Read functional annotation
# ==========================================================================

anno <- read_tsv(
  annotation_file,
  show_col_types = FALSE,
  progress = FALSE
)

required_columns <- c(
  "gene_id",
  "eggnog_cog_category",
  "signalp_is_sp",
  "dbcan_high_confidence"
)

missing_columns <- setdiff(required_columns, names(anno))

if (length(missing_columns) > 0) {
  stop(
    "Missing annotation column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

if (anyDuplicated(anno$gene_id)) {
  stop("Publication annotation contains duplicated gene IDs.")
}

# ==========================================================================
# Join
# ==========================================================================

x <- spatial |>
  left_join(
    anno |>
      select(all_of(required_columns)),
    by = "gene_id"
  )

if (nrow(x) != nrow(spatial)) {
  stop("Join changed the number of classified genes.")
}

missing_annotation <- sum(
  is.na(x$signalp_is_sp) &
    is.na(x$dbcan_high_confidence) &
    is.na(x$eggnog_cog_category)
)

if (missing_annotation > 0) {
  stop(
    "Missing functional annotation for ",
    missing_annotation,
    " classified genes."
  )
}

x <- x |>
  mutate(
    is_secreted = truthy(signalp_is_sp),
    is_cazyme = truthy(dbcan_high_confidence),
    cog_codes = map(
      as.character(eggnog_cog_category),
      function(z) {
        if (is.na(z) || z %in% c("", "-")) {
          return(character())
        }

        codes <- str_extract_all(
          str_to_upper(z),
          "[A-Z]"
        )[[1]]

        sort(
          unique(
            codes[codes %in% valid_codes]
          )
        )
      }
    )
  )

# ==========================================================================
# Descriptive functional summary
# ==========================================================================

functional_summary <- x |>
  group_by(response_group) |>
  summarise(
    genes = n(),
    secreted_genes = sum(is_secreted),
    secreted_percent = 100 * mean(is_secreted),
    high_confidence_cazyme_genes = sum(is_cazyme),
    cazyme_percent = 100 * mean(is_cazyme),
    .groups = "drop"
  ) |>
  arrange(match(response_group, expected_groups))

write_tsv(
  functional_summary,
  summary_out
)

# ==========================================================================
# Secretome / CAZyme enrichment
# ==========================================================================

test_binary_feature <- function(feature, feature_name) {

  map_dfr(
    expected_groups,
    function(group_name) {

      in_group <- x$response_group == group_name
      positive <- x[[feature]]

      contingency <- matrix(
        c(
          sum(in_group & positive),
          sum(in_group & !positive),
          sum(!in_group & positive),
          sum(!in_group & !positive)
        ),
        nrow = 2,
        byrow = TRUE
      )

      ft <- fisher.test(contingency)

      tibble(
        response_group = group_name,
        feature = feature_name,
        group_positive = sum(in_group & positive),
        group_total = sum(in_group),
        group_percent = 100 * mean(positive[in_group]),
        other_positive = sum(!in_group & positive),
        other_total = sum(!in_group),
        other_percent = 100 * mean(positive[!in_group]),
        odds_ratio = unname(ft$estimate),
        pvalue = ft$p.value
      )
    }
  )
}

binary_results <- bind_rows(
  test_binary_feature(
    "is_secreted",
    "Secreted protein"
  ),
  test_binary_feature(
    "is_cazyme",
    "High-confidence CAZyme"
  )
) |>
  group_by(feature) |>
  mutate(
    padj = p.adjust(pvalue, method = "BH")
  ) |>
  ungroup() |>
  mutate(
    significant = padj < alpha,
    enrichment_direction = case_when(
      !significant ~ "Not significant",
      odds_ratio > 1 ~ "Enriched",
      odds_ratio < 1 ~ "Depleted",
      TRUE ~ "No difference"
    )
  ) |>
  arrange(feature, padj)

write_tsv(
  binary_results,
  binary_out
)

# ==========================================================================
# COG enrichment
# ==========================================================================

cog_long <- x |>
  select(
    gene_id,
    response_group,
    cog_codes
  ) |>
  unnest(cog_codes) |>
  distinct()

cog_results <- map_dfr(
  expected_groups,
  function(group_name) {

    map_dfr(
      valid_codes,
      function(code) {

        genes_with_code <- cog_long |>
          filter(cog_codes == code) |>
          pull(gene_id)

        in_group <- x$response_group == group_name
        has_code <- x$gene_id %in% genes_with_code

        a <- sum(in_group & has_code)
        b <- sum(in_group & !has_code)
        c <- sum(!in_group & has_code)
        d <- sum(!in_group & !has_code)

        if ((a + c) == 0) {
          return(NULL)
        }

        ft <- fisher.test(
          matrix(
            c(a, b, c, d),
            nrow = 2,
            byrow = TRUE
          )
        )

        tibble(
          response_group = group_name,
          cog = code,
          group_hits = a,
          group_total = sum(in_group),
          group_percent = 100 * a / sum(in_group),
          other_hits = c,
          other_total = sum(!in_group),
          other_percent = 100 * c / sum(!in_group),
          odds_ratio = unname(ft$estimate),
          pvalue = ft$p.value
        )
      }
    )
  }
) |>
  group_by(response_group) |>
  mutate(
    padj = p.adjust(pvalue, method = "BH")
  ) |>
  ungroup() |>
  mutate(
    significant = padj < alpha,
    enrichment_direction = case_when(
      !significant ~ "Not significant",
      odds_ratio > 1 ~ "Enriched",
      odds_ratio < 1 ~ "Depleted",
      TRUE ~ "No difference"
    )
  ) |>
  arrange(response_group, padj)

write_tsv(
  cog_results,
  cog_complete_out
)

write_tsv(
  cog_results |>
    filter(significant),
  cog_significant_out
)

# ==========================================================================
# Run information
# ==========================================================================

run_info <- tibble(
  field = c(
    "analysis",
    "classification_file",
    "annotation_file",
    "classified_gene_count",
    "alpha",
    "multiple_testing_binary",
    "multiple_testing_cog",
    "comparison"
  ),
  value = c(
    "Ophiostoma spatial functional enrichment",
    classification_file,
    annotation_file,
    as.character(nrow(x)),
    as.character(alpha),
    "BH correction within feature across four spatial programs",
    "BH correction within spatial program across COG categories",
    "each spatial response program versus all other classified genes"
  )
)

write_tsv(
  run_info,
  run_info_out
)

# ==========================================================================
# Console report
# ==========================================================================

cat("\n============================================================\n")
cat("SPATIAL PROGRAM FUNCTIONAL SUMMARY\n")
cat("============================================================\n")

print(
  functional_summary,
  n = Inf,
  width = Inf
)

cat("\n============================================================\n")
cat("SIGNIFICANT SECRETOME / CAZYME RESULTS\n")
cat("============================================================\n")

print(
  binary_results |>
    filter(significant) |>
    select(
      response_group,
      feature,
      group_positive,
      group_total,
      group_percent,
      other_percent,
      odds_ratio,
      padj,
      enrichment_direction
    ),
  n = Inf,
  width = Inf
)

cat("\n============================================================\n")
cat("SIGNIFICANT COG RESULTS\n")
cat("============================================================\n")

print(
  cog_results |>
    filter(significant) |>
    select(
      response_group,
      cog,
      group_hits,
      group_total,
      group_percent,
      other_percent,
      odds_ratio,
      padj,
      enrichment_direction
    ),
  n = Inf,
  width = Inf
)

cat("\nPASS: spatial functional enrichment analysis complete.\n")
cat("Outputs:\n")
cat("  ", summary_out, "\n")
cat("  ", binary_out, "\n")
cat("  ", cog_complete_out, "\n")
cat("  ", cog_significant_out, "\n")
cat("  ", run_info_out, "\n")
