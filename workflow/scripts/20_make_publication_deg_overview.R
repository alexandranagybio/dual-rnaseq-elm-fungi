#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

# -------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------

fusarium_file <- paste0(
  "results/fusarium/deseq2_results/",
  "fusarium_interaction_vs_self_raw.tsv"
)

ophiostoma_root <- "results/ophiostoma/deseq2_results"

output_dir <- "figures/publication"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_pdf <- file.path(
  output_dir,
  "figure_3c_extent_direction_differential_expression.pdf"
)

output_png <- file.path(
  output_dir,
  "figure_3c_extent_direction_differential_expression.png"
)

results_dir <- file.path(
  "results",
  "publication",
  "figure3",
  "extent_direction"
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_tsv <- file.path(
  results_dir,
  "figure_3c_extent_direction_differential_expression_data.tsv"
)

# remains-inspired organism colours
fusarium_colour <- "#E79A4B"
ophiostoma_colour <- "#8B6A9E"

padj_cutoff <- 0.05

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

validate_deseq_table <- function(x, source_name) {
  required <- c("log2FoldChange", "padj")
  missing <- setdiff(required, names(x))

  if (length(missing) > 0) {
    stop(
      "Missing required column(s) in ", source_name, ": ",
      paste(missing, collapse = ", ")
    )
  }

  x
}

summarise_contrast <- function(
    x,
    organism,
    contrast,
    display_label,
    source_file
) {
  x <- validate_deseq_table(x, source_file)

  # Genes with non-NA adjusted P-values form the denominator used here.
  tested <- sum(!is.na(x$padj))

  significant_up <- sum(
    !is.na(x$padj) &
      x$padj < padj_cutoff &
      !is.na(x$log2FoldChange) &
      x$log2FoldChange > 0
  )

  significant_down <- sum(
    !is.na(x$padj) &
      x$padj < padj_cutoff &
      !is.na(x$log2FoldChange) &
      x$log2FoldChange < 0
  )

  if (tested == 0) {
    stop("No genes with non-NA padj values in: ", source_file)
  }

  tibble(
    organism = organism,
    contrast = contrast,
    display_label = display_label,
    tested_genes = tested,
    significant_down = significant_down,
    significant_up = significant_up,
    significant_total = significant_down + significant_up,
    percent_down = 100 * significant_down / tested,
    percent_up = 100 * significant_up / tested,
    percent_total = 100 * (significant_down + significant_up) / tested,
    source_file = source_file
  )
}

read_deseq_table <- function(path) {
  read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

# Locate one Ophiostoma contrast table using filename patterns and verify
# that it contains the required DESeq2 columns.
find_ophiostoma_table <- function(include_patterns, exclude_patterns = NULL) {
  candidates <- list.files(
    ophiostoma_root,
    pattern = "\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )

  for (pattern in include_patterns) {
    candidates <- candidates[
      str_detect(
        str_to_lower(basename(candidates)),
        str_to_lower(pattern)
      )
    ]
  }

  if (!is.null(exclude_patterns)) {
    for (pattern in exclude_patterns) {
      candidates <- candidates[
        !str_detect(
          str_to_lower(basename(candidates)),
          str_to_lower(pattern)
        )
      ]
    }
  }

  # Ignore obvious summaries and diagnostics.
  candidates <- candidates[
    !str_detect(
      str_to_lower(candidates),
      "summary|validation|diagnostic|run_info|annotation"
    )
  ]

  valid <- candidates[vapply(
    candidates,
    function(path) {
      header <- tryCatch(
        names(read_tsv(
          path,
          n_max = 1,
          show_col_types = FALSE,
          progress = FALSE
        )),
        error = function(e) character()
      )

      all(c("log2FoldChange", "padj") %in% header)
    },
    logical(1)
  )]

  if (length(valid) != 1) {
    stop(
      "\nCould not uniquely identify an Ophiostoma table.\n",
      "Patterns: ", paste(include_patterns, collapse = ", "), "\n",
      "Matching valid files:\n",
      paste0("  - ", valid, collapse = "\n")
    )
  }

  valid
}

# -------------------------------------------------------------------------
# Find input tables
# -------------------------------------------------------------------------

if (!file.exists(fusarium_file)) {
  stop("Missing Fusarium DESeq2 table: ", fusarium_file)
}

ophiostoma_interaction_vs_self <- file.path(
  ophiostoma_root,
  "tables",
  "interaction_vs_self_all_genes.tsv"
)

ophiostoma_onu_vs_self <- file.path(
  ophiostoma_root,
  "tables",
  "onu_vs_self_all_genes.tsv"
)

ophiostoma_interaction_vs_onu <- file.path(
  ophiostoma_root,
  "tables",
  "interaction_vs_onu_all_genes.tsv"
)

ophiostoma_files <- c(
  ophiostoma_interaction_vs_self,
  ophiostoma_onu_vs_self,
  ophiostoma_interaction_vs_onu
)

missing_files <- ophiostoma_files[!file.exists(ophiostoma_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing Ophiostoma DESeq2 table(s):\n",
    paste0("  - ", missing_files, collapse = "\n")
  )
}

# -------------------------------------------------------------------------
# Calculate counts and percentages
# -------------------------------------------------------------------------

summary_data <- bind_rows(
  summarise_contrast(
    read_deseq_table(fusarium_file),
    organism = "Fusarium",
    contrast = "interaction_vs_self",
    display_label = expression(
      italic("Fusarium") * ": interaction vs control"
    ),
    source_file = fusarium_file
  ),
  summarise_contrast(
    read_deseq_table(ophiostoma_interaction_vs_onu),
    organism = "Ophiostoma",
    contrast = "interaction_vs_onu",
    display_label = expression(
      italic("Ophiostoma") * ": interaction vs control"
    ),
    source_file = ophiostoma_interaction_vs_onu
  ),
  summarise_contrast(
    read_deseq_table(ophiostoma_onu_vs_self),
    organism = "Ophiostoma",
    contrast = "onu_vs_self",
    display_label = expression(
      italic("Ophiostoma") * ": control vs distal"
    ),
    source_file = ophiostoma_onu_vs_self
  ),
  summarise_contrast(
    read_deseq_table(ophiostoma_interaction_vs_self),
    organism = "Ophiostoma",
    contrast = "interaction_vs_self",
    display_label = expression(
      italic("Ophiostoma") * ": interaction vs distal"
    ),
    source_file = ophiostoma_interaction_vs_self
  )
)

write_tsv(
  summary_data |>
    select(-display_label),
  output_tsv
)

print(
  summary_data |>
    select(
      organism,
      contrast,
      tested_genes,
      significant_down,
      significant_up,
      significant_total,
      percent_down,
      percent_up,
      percent_total
    )
)

# -------------------------------------------------------------------------
# Prepare plotting data
# -------------------------------------------------------------------------

# Figure 3C compares the biologically equivalent confrontation contrast
# for each organism. All four contrasts remain available in the output TSV.
plot_data <- summary_data |>
  filter(
    (
      organism == "Fusarium" &
        contrast == "interaction_vs_self"
    ) |
      (
        organism == "Ophiostoma" &
          contrast == "interaction_vs_onu"
      )
  ) |>
  mutate(
    plot_label = case_when(
      organism == "Fusarium" ~ "Fusarium cf. salinense",
      organism == "Ophiostoma" ~ "Ophiostoma novo-ulmi"
    ),
    plot_label = factor(
      plot_label,
      levels = c(
        "Ophiostoma novo-ulmi",
        "Fusarium cf. salinense"
      )
    )
  ) |>
  select(
    organism,
    contrast,
    plot_label,
    significant_down,
    significant_up
  ) |>
  pivot_longer(
    cols = c(significant_down, significant_up),
    names_to = "direction",
    values_to = "count"
  ) |>
  mutate(
    signed_count = if_else(
      direction == "significant_down",
      -count,
      count
    )
  )

# Symmetric count axis with sufficient room for the end labels.
largest_extent <- max(abs(plot_data$signed_count))
axis_limit <- ceiling((largest_extent + 300) / 500) * 500
axis_breaks <- seq(-axis_limit, axis_limit, by = 1000)

# -------------------------------------------------------------------------
# Plot
# -------------------------------------------------------------------------

p <- ggplot(
  plot_data,
  aes(
    x = signed_count,
    y = plot_label,
    fill = organism
  )
) +
  geom_col(
    width = 0.64,
    alpha = 0.70
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.42,
    colour = "grey72"
  ) +
  geom_text(
    aes(
      label = scales::comma(count),
      hjust = if_else(signed_count < 0, 1.08, -0.08)
    ),
    size = 3.5,
    family = "Nimbus Sans",
    colour = "grey15"
  ) +
  scale_fill_manual(
    values = c(
      Fusarium = fusarium_colour,
      Ophiostoma = ophiostoma_colour
    ),
    guide = "none"
  ) +
  scale_y_discrete(
    labels = expression(
      italic("Ophiostoma novo-ulmi"),
      italic("Fusarium cf. salinense")
    )
  ) +
  scale_x_continuous(
    limits = c(-axis_limit, axis_limit),
    breaks = axis_breaks,
    labels = function(x) scales::comma(abs(x)),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Significantly differentially expressed genes",
    y = NULL,
    tag = "c"
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  theme_classic(base_size = 10, base_family = "Nimbus Sans") +
  theme(
    axis.text.x = element_text(
      colour = "grey20",
      size = 9
    ),
    axis.text.y = element_text(
      colour = "grey15",
      size = 10,
      margin = margin(r = 10)
    ),
    axis.title.x = element_text(
      size = 10,
      margin = margin(t = 8)
    ),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    plot.tag = element_text(
      face = "bold",
      size = 10,
      colour = "grey10"
    ),
    plot.tag.position = c(0, 1),
    plot.margin = margin(13, 28, 8, 14)
  )

ggsave(
  output_pdf,
  plot = p,
  width = 170,
  height = 52.6,
  units = "mm",
  device = cairo_pdf
)

ggsave(
  output_png,
  plot = p,
  width = 170,
  height = 52.6,
  units = "mm",
  dpi = 600,
  bg = "white"
)

message("\nCreated:")
message("  ", output_pdf)
message("  ", output_png)
message("  ", output_tsv)
