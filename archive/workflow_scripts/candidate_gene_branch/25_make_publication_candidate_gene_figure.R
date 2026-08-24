#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(patchwork)
})

# ==========================================================================
# Configuration
# ==========================================================================

n_per_direction <- 8
alpha_threshold <- 0.05
lfc_threshold <- 1

fusarium_colour <- "#D9792B"
ophiostoma_colour <- "#76519D"

output_dir <- "figures/figure4_candidate_genes"
table_dir <- file.path(
  "results/publication",
  "figure4_candidate_genes"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

fusarium_input <- paste0(
  "results/fusarium/publication_annotation/tables/",
  "fusarium_publication_annotation_significant.tsv"
)

ophiostoma_input <- paste0(
  "results/ophiostoma/publication_annotation/tables/",
  "interaction_vs_self_publication_annotation.tsv"
)

ophiostoma_shrunk_input <- paste0(
  "results/ophiostoma/deseq2_results/tables/",
  "interaction_vs_self_lfcshrunk.tsv"
)

# ==========================================================================
# Helpers
# ==========================================================================

clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- str_squish(x)

  missing_values <- c(
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

  x[x %in% missing_values] <- ""
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

  output <- rep(
    "",
    length(fields[[1]])
  )

  for (field in fields) {
    use <- output == "" & field != ""
    output[use] <- field[use]
  }

  output
}

sentence_case_safe <- function(x) {
  x <- clean_text(x)

  ifelse(
    x == "",
    "",
    paste0(
      str_to_upper(str_sub(x, 1, 1)),
      str_sub(x, 2)
    )
  )
}

shorten_label <- function(x, width = 58) {
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

is_generic_annotation <- function(x) {
  value <- str_to_lower(
    clean_text(x)
  )

  value == "" |
    str_detect(
      value,
      paste(
        c(
          "^hypothetical protein",
          "^uncharacteri[sz]ed protein",
          "^predicted protein",
          "^protein of unknown function",
          "^unknown protein",
          "^conserved hypothetical protein",
          "^domain of unknown function",
          "^integral membrane protein$",
          "^membrane protein$",
          "^zinc finger$",
          "^gata zinc finger$",
          "^fungal specific transcription factor domain$",
          "^protein kinase domain$",
          "^cazyme family$",
          "^acting on the ch-oh group of donors",
          "^oxidoreductase activity$",
          "^hydrolase activity$",
          "^binding$",
          "^catalytic activity$",
          "^metabolic process$",
          "^cellular process$"
        ),
        collapse = "|"
      )
    )
}

prepare_labels <- function(
    preferred_name,
    description,
    pfam,
    gene_id
) {
  preferred_name <- clean_text(preferred_name)
  description <- clean_text(description)
  pfam <- clean_text(pfam)

  preferred_name[
    is_generic_annotation(preferred_name)
  ] <- ""

  description[
    is_generic_annotation(description)
  ] <- ""

  pfam[
    is_generic_annotation(pfam)
  ] <- ""

  label <- first_nonempty(
    preferred_name,
    description,
    pfam
  )

  label <- sentence_case_safe(label)
  label <- shorten_label(label)

  label
}

select_candidates <- function(
    dat,
    n_each = 8
) {
  dat %>%
    filter(
      informative,
      !is.na(plot_lfc),
      !is.na(padj),
      padj < alpha_threshold,
      abs(plot_lfc) > lfc_threshold
    ) %>%
    arrange(
      direction,
      desc(abs(plot_lfc)),
      padj,
      desc(base_mean)
    ) %>%
    group_by(
      organism,
      direction,
      display_label
    ) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    group_by(
      organism,
      direction
    ) %>%
    slice_head(n = n_each) %>%
    ungroup()
}

make_panel <- function(
    dat,
    species_colour,
    panel_letter,
    species_label
) {
  dat <- dat %>%
    arrange(plot_lfc) %>%
    mutate(
      row_label = paste0(
        display_label,
        "  [",
        gene_id,
        "]"
      ),
      row_label = factor(
        row_label,
        levels = row_label
      )
    )

  x_limit <- max(
    abs(dat$plot_lfc),
    na.rm = TRUE
  )

  x_limit <- ceiling(
    x_limit * 1.12
  )

  ggplot(
    dat,
    aes(
      x = plot_lfc,
      y = row_label
    )
  ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.35,
      colour = "grey70"
    ) +
    geom_segment(
      aes(
        x = 0,
        xend = plot_lfc,
        yend = row_label
      ),
      linewidth = 0.75,
      colour = species_colour,
      alpha = 0.65
    ) +
    geom_point(
      size = 3.1,
      colour = species_colour
    ) +
    geom_text(
      aes(
        label = sprintf(
          "%.2f",
          plot_lfc
        )
      ),
      hjust = ifelse(
        dat$plot_lfc > 0,
        -0.28,
        1.28
      ),
      size = 3.1,
      colour = "grey25"
    ) +
    scale_x_continuous(
      limits = c(
        -x_limit,
        x_limit
      ),
      breaks = scales::pretty_breaks(n = 5),
      expand = expansion(
        mult = c(0.04, 0.04)
      )
    ) +
    labs(
      title = paste0(
        panel_letter,
        "   ",
        species_label
      ),
      x = expression(log[2]~fold~change),
      y = NULL
    ) +
    theme_classic(
      base_size = 10.5
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 11.5,
        margin = margin(
          b = 8
        )
      ),
      axis.text.y = element_text(
        size = 8.4,
        colour = "grey15"
      ),
      axis.text.x = element_text(
        colour = "grey25"
      ),
      axis.title.x = element_text(
        margin = margin(
          t = 7
        )
      ),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(
        8,
        12,
        8,
        8
      )
    )
}

# ==========================================================================
# Fusarium
# ==========================================================================

fusarium_raw <- read_tsv(
  fusarium_input,
  show_col_types = FALSE,
  progress = FALSE
)

fusarium <- fusarium_raw %>%
  transmute(
    organism = "Fusarium",
    gene_id,
    base_mean = raw_baseMean,
    raw_lfc = raw_log2FoldChange,
    plot_lfc = shrunk_log2FoldChange,
    padj = raw_padj,

    preferred_name = clean_text(
      eggnog_preferred_name
    ),

    description = clean_text(
      eggnog_description
    ),

    pfam = clean_text(
      eggnog_pfams
    ),

    display_label = prepare_labels(
      preferred_name,
      description,
      pfam,
      gene_id
    ),

    informative = (
      display_label != ""
    ),

    direction = if_else(
      plot_lfc > 0,
      "Upregulated",
      "Downregulated"
    ),

    eggnog_annotated = as_logical_safe(
      eggnog_annotated
    ),

    signal_peptide = as_logical_safe(
      signalp_confident_positive
    ),

    cazyme = as_logical_safe(
      dbcan_any_hit
    )
  )

fusarium_selected <- select_candidates(
  fusarium,
  n_per_direction
)

# ==========================================================================
# Ophiostoma
# ==========================================================================

ophiostoma_raw <- read_tsv(
  ophiostoma_input,
  show_col_types = FALSE,
  progress = FALSE
)

ophiostoma_shrunk <- read_tsv(
  ophiostoma_shrunk_input,
  show_col_types = FALSE,
  progress = FALSE
) %>%
  select(
    gene_id,
    shrunk_log2FoldChange,
    shrunk_lfcSE
  )

if (anyDuplicated(ophiostoma_shrunk$gene_id) > 0) {
  stop(
    "Duplicated gene IDs in Ophiostoma shrinkage table."
  )
}

ophiostoma_joined <- ophiostoma_raw %>%
  left_join(
    ophiostoma_shrunk,
    by = "gene_id"
  )

if (nrow(ophiostoma_joined) != nrow(ophiostoma_raw)) {
  stop(
    "Ophiostoma row count changed after joining shrinkage estimates."
  )
}

if (
  any(
    is.na(
      ophiostoma_joined$shrunk_log2FoldChange
    )
  )
) {
  stop(
    "Missing Ophiostoma shrunken LFC values after join."
  )
}

ophiostoma <- ophiostoma_joined %>%
  filter(
    as_logical_safe(
      de_significant_padj_lt_0.05
    )
  ) %>%
  transmute(
    organism = "Ophiostoma",
    gene_id,
    base_mean = baseMean,
    raw_lfc = log2FoldChange,
    plot_lfc = shrunk_log2FoldChange,
    padj,

    preferred_name = clean_text(
      eggnog_preferred_name
    ),

    description = clean_text(
      eggnog_description
    ),

    pfam = clean_text(
      eggnog_pfams
    ),

    display_label = prepare_labels(
      preferred_name,
      description,
      pfam,
      gene_id
    ),

    informative = (
      display_label != ""
    ),

    direction = if_else(
      plot_lfc > 0,
      "Upregulated",
      "Downregulated"
    ),

    eggnog_annotated = as_logical_safe(
      eggnog_annotated
    ),

    signal_peptide = as_logical_safe(
      signalp_is_sp
    ),

    cazyme = as_logical_safe(
      dbcan_any_hit
    )
  )

ophiostoma_selected <- select_candidates(
  ophiostoma,
  n_per_direction
)

# ==========================================================================
# Validation
# ==========================================================================

expected_per_species <- (
  n_per_direction * 2
)

if (
  nrow(fusarium_selected)
  != expected_per_species
) {
  stop(
    "Fusarium selection produced ",
    nrow(fusarium_selected),
    " candidates; expected ",
    expected_per_species
  )
}

if (
  nrow(ophiostoma_selected)
  != expected_per_species
) {
  stop(
    "Ophiostoma selection produced ",
    nrow(ophiostoma_selected),
    " candidates; expected ",
    expected_per_species
  )
}

selection_summary <- bind_rows(
  fusarium_selected,
  ophiostoma_selected
) %>%
  count(
    organism,
    direction,
    name = "selected_genes"
  )

# ==========================================================================
# Write selected candidate tables
# ==========================================================================

write_tsv(
  fusarium_selected %>%
    arrange(
      direction,
      desc(abs(plot_lfc))
    ),
  file.path(
    table_dir,
    "fusarium_selected_candidate_genes.tsv"
  )
)

write_tsv(
  ophiostoma_selected %>%
    arrange(
      direction,
      desc(abs(plot_lfc))
    ),
  file.path(
    table_dir,
    "ophiostoma_selected_candidate_genes.tsv"
  )
)

write_tsv(
  bind_rows(
    fusarium_selected,
    ophiostoma_selected
  ) %>%
    arrange(
      organism,
      direction,
      desc(abs(plot_lfc))
    ),
  file.path(
    table_dir,
    "figure4_selected_candidate_genes.tsv"
  )
)

write_tsv(
  selection_summary,
  file.path(
    table_dir,
    "figure4_selection_summary.tsv"
  )
)

# ==========================================================================
# Plot
# ==========================================================================

panel_a <- make_panel(
  fusarium_selected,
  fusarium_colour,
  "A",
  expression(
    italic("Fusarium cf. salinense")
  )
)

panel_b <- make_panel(
  ophiostoma_selected,
  ophiostoma_colour,
  "B",
  expression(
    italic("Ophiostoma novo-ulmi")
  )
)

figure_4 <- panel_a + panel_b +
  plot_layout(
    widths = c(1, 1)
  )

ggsave(
  filename = file.path(
    output_dir,
    "Figure4_candidate_genes_draft.pdf"
  ),
  plot = figure_4,
  width = 13.5,
  height = 7.7,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure4_candidate_genes_draft.png"
  ),
  plot = figure_4,
  width = 13.5,
  height = 7.7,
  units = "in",
  dpi = 400,
  bg = "white"
)

cat("\n")
cat("============================================================\n")
cat("FIGURE 4 CANDIDATE-GENE DRAFT COMPLETE\n")
cat("============================================================\n")

cat("\nSelection summary:\n")
print(
  selection_summary,
  n = Inf
)

cat("\nFusarium candidates:\n")
fusarium_selected %>%
  select(
    gene_id,
    direction,
    plot_lfc,
    padj,
    display_label
  ) %>%
  arrange(
    direction,
    desc(abs(plot_lfc))
  ) %>%
  print(
    n = Inf,
    width = Inf
  )

cat("\nOphiostoma candidates:\n")
ophiostoma_selected %>%
  select(
    gene_id,
    direction,
    plot_lfc,
    padj,
    display_label
  ) %>%
  arrange(
    direction,
    desc(abs(plot_lfc))
  ) %>%
  print(
    n = Inf,
    width = Inf
  )

cat("\nOutputs:\n")
cat(
  normalizePath(output_dir),
  "\n"
)

cat(
  normalizePath(table_dir),
  "\n"
)
