#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

out_dir <- file.path(
  "results/publication",
  "figure4_candidate_genes",
  "candidate_dossier"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

fusarium_path <- paste0(
  "results/fusarium/publication_annotation/tables/",
  "fusarium_publication_annotation_significant.tsv"
)

ophiostoma_annotation_path <- paste0(
  "results/ophiostoma/publication_annotation/tables/",
  "interaction_vs_self_publication_annotation.tsv"
)

ophiostoma_shrunk_path <- paste0(
  "results/ophiostoma/deseq2_results/tables/",
  "interaction_vs_self_lfcshrunk.tsv"
)

selected_path <- paste0(
  "results/publication/figure4_candidate_genes/",
  "figure4_selected_candidate_genes.tsv"
)

alpha_threshold <- 0.05
lfc_threshold <- 1
pool_per_direction <- 30

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- str_squish(x)

  x[
    x %in% c(
      "",
      "-",
      "--",
      ".",
      "NA",
      "N/A",
      "na",
      "None",
      "none",
      "unknown",
      "Unknown"
    )
  ] <- ""

  x
}

as_logical_safe <- function(x) {
  if (is.logical(x)) {
    x[is.na(x)] <- FALSE
    return(x)
  }

  out <- str_to_lower(as.character(x)) %in%
    c("true", "t", "1", "yes", "y")

  out[is.na(out)] <- FALSE
  out
}

first_nonempty <- function(...) {
  fields <- lapply(
    list(...),
    clean_text
  )

  result <- rep(
    "",
    length(fields[[1]])
  )

  for (field in fields) {
    use <- result == "" & field != ""
    result[use] <- field[use]
  }

  result
}

truncate_text <- function(x, width = 160) {
  x <- clean_text(x)

  ifelse(
    nchar(x) > width,
    paste0(
      str_sub(x, 1, width - 1),
      "…"
    ),
    x
  )
}

is_generic <- function(x) {
  x <- str_to_lower(
    clean_text(x)
  )

  x == "" |
    str_detect(
      x,
      paste(
        c(
          "^hypothetical protein",
          "^uncharacteri[sz]ed",
          "^predicted protein",
          "^protein of unknown function",
          "^unknown protein",
          "^conserved hypothetical",
          "^integral membrane protein$",
          "^membrane protein$",
          "^domain of unknown function",
          "^fungal_trans$",
          "^bzip_2$",
          "^act domain$",
          "^aaa domain$",
          "^fad binding domain$",
          "^bcs1 n terminal$"
        ),
        collapse = "|"
      )
    )
}

assign_theme <- function(
    preferred_name,
    description,
    pfam,
    cog,
    kegg,
    go,
    cazy,
    dbcan
) {
  text <- str_to_lower(
    paste(
      clean_text(preferred_name),
      clean_text(description),
      clean_text(pfam),
      clean_text(cog),
      clean_text(kegg),
      clean_text(go),
      clean_text(cazy),
      clean_text(dbcan),
      collapse = " | "
    )
  )

  case_when(
    str_detect(
      text,
      "transcription factor|zinc cluster|bzip|gata|receptor|signalling|signaling|protein kinase"
    ) ~ "Regulation/signalling",

    str_detect(
      text,
      "transporter|permease|abc-|major facilitator|mfs|amino acid transport|nitrate|sulfonate|bicarbonate"
    ) ~ "Transport",

    str_detect(
      text,
      "oxidoreductase|dehydrogenase|reductase|oxygenase|p450|cytochrome|aldo/keto|quinone|redox"
    ) ~ "Redox/detoxification",

    str_detect(
      text,
      "cellulose|glucan|glycoside hydrolase|cazy|carbohydrate-active|cell wall|acid phosphatase|bgl"
    ) ~ "Cell wall/extracellular",

    str_detect(
      text,
      "acetyl-coa|amino transferase|aminotransferase|5-oxoprolinase|metabolism|biosynthesis|catabolism|amp-binding|fatty acid"
    ) ~ "Central/specialized metabolism",

    str_detect(
      text,
      "dna repair|base excision|ankyrin|ubiquitin|proteostasis|aaa"
    ) ~ "Stress/proteostasis",

    TRUE ~ "Unresolved"
  )
}

annotation_quality <- function(
    preferred_name,
    description,
    pfam
) {
  preferred_name <- clean_text(preferred_name)
  description <- clean_text(description)
  pfam <- clean_text(pfam)

  case_when(
    preferred_name != "" &
      !is_generic(preferred_name) ~ "Named gene",

    description != "" &
      !is_generic(description) ~ "Specific description",

    pfam != "" &
      !is_generic(pfam) ~ "Domain-supported",

    TRUE ~ "Weak/generic"
  )
}

# -------------------------------------------------------------------------
# Fusarium
# -------------------------------------------------------------------------

fusarium <- read_tsv(
  fusarium_path,
  show_col_types = FALSE,
  progress = FALSE
) %>%
  transmute(
    organism = "Fusarium",
    gene_id,
    base_mean = raw_baseMean,
    raw_lfc = raw_log2FoldChange,
    shrunk_lfc = shrunk_log2FoldChange,
    padj = raw_padj,
    direction = if_else(
      shrunk_lfc > 0,
      "Upregulated",
      "Downregulated"
    ),

    eggnog_annotated = as_logical_safe(
      eggnog_annotated
    ),

    preferred_name = clean_text(
      eggnog_preferred_name
    ),

    description = clean_text(
      eggnog_description
    ),

    seed_ortholog = clean_text(
      eggnog_seed_ortholog
    ),

    ortholog_groups = clean_text(
      eggnog_eggnog_ogs
    ),

    max_annotation_level = clean_text(
      eggnog_max_annot_lvl
    ),

    cog_category = clean_text(
      eggnog_cog_category
    ),

    pfam = clean_text(
      eggnog_pfams
    ),

    eggnog_go = clean_text(
      eggnog_gos
    ),

    eggnog_ec = clean_text(
      eggnog_ec
    ),

    eggnog_kegg_ko = clean_text(
      eggnog_kegg_ko
    ),

    eggnog_kegg_pathway = clean_text(
      eggnog_kegg_pathway
    ),

    eggnog_cazy = clean_text(
      eggnog_cazy
    ),

    signal_peptide = as_logical_safe(
      signalp_confident_positive
    ),

    dbcan_any_hit = as_logical_safe(
      dbcan_any_hit
    ),

    dbcan_high_confidence = as_logical_safe(
      dbcan_high_confidence
    ),

    dbcan_family = clean_text(
      dbcan_recommended
    ),

    dbcan_substrate = clean_text(
      dbcan_substrate
    )
  )

# -------------------------------------------------------------------------
# Ophiostoma
# -------------------------------------------------------------------------

ophiostoma_annotation <- read_tsv(
  ophiostoma_annotation_path,
  show_col_types = FALSE,
  progress = FALSE
)

ophiostoma_shrunk <- read_tsv(
  ophiostoma_shrunk_path,
  show_col_types = FALSE,
  progress = FALSE
) %>%
  select(
    gene_id,
    shrunk_log2FoldChange,
    shrunk_lfcSE
  )

ophiostoma <- ophiostoma_annotation %>%
  left_join(
    ophiostoma_shrunk,
    by = "gene_id"
  ) %>%
  transmute(
    organism = "Ophiostoma",
    gene_id,
    base_mean = baseMean,
    raw_lfc = log2FoldChange,
    shrunk_lfc = shrunk_log2FoldChange,
    padj,
    direction = if_else(
      shrunk_lfc > 0,
      "Upregulated",
      "Downregulated"
    ),

    eggnog_annotated = as_logical_safe(
      eggnog_annotated
    ),

    preferred_name = clean_text(
      eggnog_preferred_name
    ),

    description = clean_text(
      eggnog_description
    ),

    seed_ortholog = clean_text(
      eggnog_seed_ortholog
    ),

    ortholog_groups = clean_text(
      eggnog_eggnog_ogs
    ),

    max_annotation_level = clean_text(
      eggnog_max_annot_lvl
    ),

    cog_category = clean_text(
      eggnog_cog_category
    ),

    pfam = clean_text(
      eggnog_pfams
    ),

    eggnog_go = clean_text(
      eggnog_gos
    ),

    eggnog_ec = clean_text(
      eggnog_ec
    ),

    eggnog_kegg_ko = clean_text(
      eggnog_kegg_ko
    ),

    eggnog_kegg_pathway = clean_text(
      eggnog_kegg_pathway
    ),

    eggnog_cazy = clean_text(
      eggnog_cazy
    ),

    jgi_go_names = clean_text(
      go_names
    ),

    jgi_kegg_definition = clean_text(
      kegg_definition
    ),

    jgi_kegg_pathway = clean_text(
      kegg_pathway
    ),

    signal_peptide = as_logical_safe(
      signalp_is_sp
    ),

    dbcan_any_hit = as_logical_safe(
      dbcan_any_hit
    ),

    dbcan_high_confidence = as_logical_safe(
      dbcan_high_confidence
    ),

    dbcan_family = clean_text(
      dbcan_recommended
    ),

    dbcan_substrate = clean_text(
      dbcan_substrate
    )
  )

# Add missing JGI columns to Fusarium so the tables bind cleanly.
fusarium <- fusarium %>%
  mutate(
    jgi_go_names = "",
    jgi_kegg_definition = "",
    jgi_kegg_pathway = ""
  )

combined <- bind_rows(
  fusarium,
  ophiostoma
) %>%
  mutate(
    display_label = first_nonempty(
      preferred_name,
      description,
      pfam,
      gene_id
    ),

    annotation_quality = annotation_quality(
      preferred_name,
      description,
      pfam
    ),

    theme = assign_theme(
      preferred_name,
      description,
      pfam,
      cog_category,
      paste(
        eggnog_kegg_ko,
        eggnog_kegg_pathway,
        jgi_kegg_definition,
        jgi_kegg_pathway
      ),
      paste(
        eggnog_go,
        jgi_go_names
      ),
      eggnog_cazy,
      paste(
        dbcan_family,
        dbcan_substrate
      )
    ),

    full_annotation = truncate_text(
      first_nonempty(
        description,
        preferred_name,
        pfam
      )
    ),

    strong = (
      !is.na(padj) &
        padj < alpha_threshold &
        !is.na(shrunk_lfc) &
        abs(shrunk_lfc) > lfc_threshold
    )
  )

# -------------------------------------------------------------------------
# Current selected candidates
# -------------------------------------------------------------------------

selected_ids <- read_tsv(
  selected_path,
  show_col_types = FALSE,
  progress = FALSE
) %>%
  select(
    organism,
    gene_id
  ) %>%
  distinct()

selected_dossier <- combined %>%
  inner_join(
    selected_ids,
    by = c(
      "organism",
      "gene_id"
    )
  ) %>%
  arrange(
    organism,
    direction,
    desc(abs(shrunk_lfc))
  )

# -------------------------------------------------------------------------
# Wider replacement pool
# -------------------------------------------------------------------------

candidate_pool <- combined %>%
  filter(
    strong,
    eggnog_annotated
  ) %>%
  arrange(
    organism,
    direction,
    factor(
      annotation_quality,
      levels = c(
        "Named gene",
        "Specific description",
        "Domain-supported",
        "Weak/generic"
      )
    ),
    desc(abs(shrunk_lfc)),
    padj,
    desc(base_mean)
  ) %>%
  group_by(
    organism,
    direction
  ) %>%
  slice_head(
    n = pool_per_direction
  ) %>%
  ungroup()

# -------------------------------------------------------------------------
# Outputs
# -------------------------------------------------------------------------

write_tsv(
  selected_dossier,
  file.path(
    out_dir,
    "current_selected_candidate_dossier.tsv"
  )
)

write_tsv(
  candidate_pool,
  file.path(
    out_dir,
    "ranked_candidate_replacement_pool.tsv"
  )
)

write_tsv(
  candidate_pool %>%
    count(
      organism,
      direction,
      annotation_quality,
      theme,
      name = "candidate_genes"
    ),
  file.path(
    out_dir,
    "candidate_pool_summary.tsv"
  )
)

cat("\n")
cat("============================================================\n")
cat("CURRENT FIGURE 4 CANDIDATE DOSSIER\n")
cat("============================================================\n")

selected_dossier %>%
  select(
    organism,
    direction,
    gene_id,
    shrunk_lfc,
    padj,
    preferred_name,
    description,
    pfam,
    cog_category,
    theme,
    annotation_quality,
    signal_peptide,
    dbcan_high_confidence
  ) %>%
  print(
    n = Inf,
    width = Inf
  )

cat("\n")
cat("============================================================\n")
cat("TOP REPLACEMENT POOL: 12 PER GROUP\n")
cat("============================================================\n")

candidate_pool %>%
  group_by(
    organism,
    direction
  ) %>%
  slice_head(n = 12) %>%
  ungroup() %>%
  select(
    organism,
    direction,
    gene_id,
    shrunk_lfc,
    preferred_name,
    description,
    pfam,
    theme,
    annotation_quality,
    signal_peptide,
    dbcan_high_confidence
  ) %>%
  print(
    n = Inf,
    width = Inf
  )

cat("\nOutputs written to:\n")
cat(normalizePath(out_dir), "\n")
