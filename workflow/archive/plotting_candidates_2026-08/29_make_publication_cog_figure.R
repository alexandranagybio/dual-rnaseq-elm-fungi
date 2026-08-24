#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

# ==========================================================================
# Configuration
# ==========================================================================

fusarium_colour <- "#D9792B"
ophiostoma_colour <- "#76519D"
neutral_colour <- "grey70"

minimum_foreground_hits <- 5
minimum_interaction_response_hits <- 20
alpha_threshold <- 0.05

enrichment_input <- file.path(
  "results",
  "publication",
  "cog_enrichment",
  "cog_enrichment_complete.tsv"
)

interaction_input <- file.path(
  "results",
  "publication",
  "cog_species_interaction",
  "species_cog_interaction_complete.tsv"
)

output_dir <- file.path(
  "figures",
  "figure4_cog_enrichment"
)

table_dir <- file.path(
  "results",
  "publication",
  "figure4_cog_final"
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

for (path in c(enrichment_input, interaction_input)) {
  if (!file.exists(path)) {
    stop("Missing required input: ", path)
  }
}

# ==========================================================================
# Read data
# ==========================================================================

enrichment <- read_tsv(
  enrichment_input,
  show_col_types = FALSE,
  progress = FALSE
)

interaction <- read_tsv(
  interaction_input,
  show_col_types = FALSE,
  progress = FALSE
)

# ==========================================================================
# Panel A: within-species enrichment
# ==========================================================================

panel_a_data <- enrichment %>%
  filter(
    cog != "S",
    foreground_hits >= minimum_foreground_hits,
    padj < alpha_threshold,
    odds_ratio > 1,
    is.finite(enrichment_ratio)
  ) %>%
  mutate(
    organism = factor(
      organism,
      levels = c(
        "Fusarium",
        "Ophiostoma"
      )
    ),

    direction = factor(
      direction,
      levels = c(
        "Induced",
        "Repressed"
      )
    ),

    category_label = paste0(
      cog,
      " · ",
      cog_name
    ),

    significance_strength = -log10(
      pmax(
        padj,
        .Machine$double.xmin
      )
    )
  )

if (nrow(panel_a_data) == 0) {
  stop("No significant within-species COG enrichments available.")
}

panel_a_order <- panel_a_data %>%
  group_by(category_label) %>%
  summarise(
    minimum_padj = min(
      padj,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(minimum_padj)
  ) %>%
  pull(category_label)

panel_a_data <- panel_a_data %>%
  mutate(
    category_label = factor(
      category_label,
      levels = panel_a_order
    )
  )

panel_a <- ggplot(
  panel_a_data,
  aes(
    x = enrichment_ratio,
    y = category_label
  )
) +
  geom_vline(
    xintercept = 1,
    linewidth = 0.35,
    colour = "grey75"
  ) +
  geom_segment(
    aes(
      x = 1,
      xend = enrichment_ratio,
      yend = category_label,
      colour = organism
    ),
    linewidth = 0.7,
    alpha = 0.55
  ) +
  geom_point(
    aes(
      size = foreground_hits,
      colour = organism,
      alpha = significance_strength
    )
  ) +
  facet_grid(
    direction ~ organism,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_colour_manual(
    values = c(
      "Fusarium" = fusarium_colour,
      "Ophiostoma" = ophiostoma_colour
    ),
    guide = "none"
  ) +
  scale_size_continuous(
    name = "Differentially expressed genes",
    range = c(3.2, 8.5)
  ) +
  scale_alpha_continuous(
    name = expression(-log[10]~adjusted~italic(P)),
    range = c(0.55, 1)
  ) +
  scale_x_continuous(
    name = "Enrichment ratio",
    expand = expansion(
      mult = c(0.03, 0.08)
    )
  ) +
  labs(
    y = NULL
  ) +
  theme_classic(
    base_size = 10.5
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 10.5
    ),
    axis.text.y = element_text(
      size = 8.2,
      colour = "grey15"
    ),
    axis.text.x = element_text(
      colour = "grey25"
    ),
    axis.title.x = element_text(
      margin = margin(t = 7)
    ),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing.x = grid::unit(
      1.25,
      "lines"
    ),
    panel.spacing.y = grid::unit(
      0.9,
      "lines"
    ),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.margin = margin(
      8,
      12,
      8,
      8
    )
  )

# ==========================================================================
# Panel B: species × COG interaction
# ==========================================================================

panel_b_data <- interaction %>%
  filter(
    cog != "S",
    interaction_padj < alpha_threshold,
    is.finite(interaction_odds_ratio),
    is.finite(interaction_ci_lower),
    is.finite(interaction_ci_upper),

    (
      fusarium_response_hits +
      ophiostoma_response_hits
    ) >= minimum_interaction_response_hits
  ) %>%
  mutate(
    response = factor(
      response,
      levels = c(
        "Induced",
        "Repressed"
      )
    ),

    category_label = paste0(
      cog,
      " · ",
      cog_name
    ),

    log2_or = log2(
      interaction_odds_ratio
    ),

    log2_ci_lower = log2(
      interaction_ci_lower
    ),

    log2_ci_upper = log2(
      interaction_ci_upper
    ),

    supported_in = case_when(
      interaction_odds_ratio < 1 ~
        "Stronger in Fusarium",

      interaction_odds_ratio > 1 ~
        "Stronger in Ophiostoma",

      TRUE ~
        "No supported difference"
    ),

    total_response_hits =
      fusarium_response_hits +
      ophiostoma_response_hits
  )

if (nrow(panel_b_data) == 0) {
  stop("No significant species × COG interactions available.")
}

panel_b_order <- panel_b_data %>%
  group_by(category_label) %>%
  summarise(
    maximum_effect = max(
      abs(log2_or),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(maximum_effect) %>%
  pull(category_label)

panel_b_data <- panel_b_data %>%
  mutate(
    category_label = factor(
      category_label,
      levels = panel_b_order
    )
  )

panel_b <- ggplot(
  panel_b_data,
  aes(
    x = log2_or,
    y = category_label
  )
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.4,
    colour = "grey65"
  ) +
  geom_errorbarh(
    aes(
      xmin = log2_ci_lower,
      xmax = log2_ci_upper,
      colour = supported_in
    ),
    height = 0,
    linewidth = 0.65,
    alpha = 0.75
  ) +
  geom_point(
    aes(
      size = total_response_hits,
      colour = supported_in
    ),
    alpha = 0.95
  ) +
  facet_wrap(
    ~ response,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_colour_manual(
    values = c(
      "Stronger in Fusarium" =
        fusarium_colour,

      "Stronger in Ophiostoma" =
        ophiostoma_colour,

      "No supported difference" =
        neutral_colour
    ),
    name = NULL
  ) +
  scale_size_continuous(
    name = "Responding genes",
    range = c(3.2, 8.5)
  ) +
  scale_x_continuous(
    name = expression(
      log[2]~
      species %*% COG~
      interaction~odds~ratio
    ),
    expand = expansion(
      mult = c(0.08, 0.08)
    )
  ) +
  labs(
    y = NULL
  ) +
  theme_classic(
    base_size = 10.5
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 10.5
    ),
    axis.text.y = element_text(
      size = 8.2,
      colour = "grey15"
    ),
    axis.text.x = element_text(
      colour = "grey25"
    ),
    axis.title.x = element_text(
      margin = margin(t = 7)
    ),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing.x = grid::unit(
      1.5,
      "lines"
    ),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.margin = margin(
      8,
      12,
      8,
      8
    )
  )

# ==========================================================================
# Combine panels
# ==========================================================================

figure_4 <- panel_a / panel_b +
  plot_layout(
    heights = c(1.08, 1)
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 14
      )
    )
  )

# ==========================================================================
# Save
# ==========================================================================

ggsave(
  filename = file.path(
    output_dir,
    "Figure4_functional_COG_programs_final.pdf"
  ),
  plot = figure_4,
  width = 13.5,
  height = 13,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure4_functional_COG_programs_final.png"
  ),
  plot = figure_4,
  width = 13.5,
  height = 13,
  units = "in",
  dpi = 400,
  bg = "white"
)

write_tsv(
  panel_a_data,
  file.path(
    table_dir,
    "figure4_panelA_within_species_enrichment.tsv"
  )
)

write_tsv(
  panel_b_data,
  file.path(
    table_dir,
    "figure4_panelB_species_cog_interactions.tsv"
  )
)

parameter_table <- tibble(
  parameter = c(
    "panel_A_threshold",
    "panel_A_minimum_foreground_hits",
    "panel_A_category_S",
    "panel_B_threshold",
    "panel_B_minimum_total_response_hits",
    "panel_B_category_S",
    "interaction_OR_interpretation"
  ),

  value = c(
    "within-species Fisher enrichment padj < 0.05 and odds ratio > 1",
    minimum_foreground_hits,
    "excluded from plotted figure",
    "species × COG binomial interaction padj < 0.05",
    minimum_interaction_response_hits,
    "excluded from plotted figure",
    "OR < 1 = stronger association in Fusarium; OR > 1 = stronger association in Ophiostoma"
  )
)

write_tsv(
  parameter_table,
  file.path(
    table_dir,
    "figure4_plot_parameters.tsv"
  )
)

cat("\n")
cat("============================================================\n")
cat("FIGURE 4 COMPLETE\n")
cat("============================================================\n")

cat("\nPanel A categories:\n")
print(
  panel_a_data %>%
    select(
      organism,
      direction,
      cog,
      cog_name,
      foreground_hits,
      enrichment_ratio,
      padj
    ),
  n = Inf,
  width = Inf
)

cat("\nPanel B significant interactions:\n")
print(
  panel_b_data %>%
    select(
      response,
      cog,
      cog_name,
      fusarium_response_hits,
      ophiostoma_response_hits,
      interaction_odds_ratio,
      interaction_ci_lower,
      interaction_ci_upper,
      interaction_padj,
      supported_in
    ),
  n = Inf,
  width = Inf
)

cat("\nFigure written to:\n")
cat(normalizePath(output_dir), "\n")

cat("\nSource-data tables written to:\n")
cat(normalizePath(table_dir), "\n")
