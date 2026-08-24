#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# ============================================================
# Configuration
# ============================================================

alpha <- 0.05

fusarium_colour <- "#D9792B"
ophiostoma_colour <- "#76519D"

enrichment_file <- file.path(
  "results",
  "publication",
  "cog_enrichment",
  "cog_enrichment_complete.tsv"
)

interaction_file <- file.path(
  "results",
  "publication",
  "cog_species_interaction",
  "species_cog_interaction_complete.tsv"
)

output_dir <- file.path(
  "figures",
  "figure4_cog_enrichment",
  "final_panels"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# Read data
# ============================================================

enrichment <- read_tsv(
  enrichment_file,
  show_col_types = FALSE
)

interaction <- read_tsv(
  interaction_file,
  show_col_types = FALSE
)

# ============================================================
# Short labels
# ============================================================

label_dictionary <- c(
  "Carbohydrate transport and metabolism" =
    "Carbohydrate metabolism",

  "Energy production and conversion" =
    "Energy production",

  "Lipid transport and metabolism" =
    "Lipid metabolism",

  "Amino-acid transport and metabolism" =
    "Amino-acid metabolism",

  "Coenzyme transport and metabolism" =
    "Coenzyme metabolism",

  "Translation and ribosome biogenesis" =
    "Translation",

  "Transcription" =
    "Transcription",

  "RNA processing and modification" =
    "RNA processing",

  "Chromatin structure and dynamics" =
    "Chromatin",

  "Post-translational modification and protein turnover" =
    "Protein turnover",

  "Intracellular trafficking and secretion" =
    "Intracellular trafficking",

  "Cytoskeleton" =
    "Cytoskeleton",

  "Cell wall, membrane and envelope biogenesis" =
    "Cell wall and membrane",

  "Cell-cycle control and chromosome partitioning" =
    "Cell cycle",

  "Secondary-metabolite biosynthesis and transport" =
    "Secondary metabolism",

  "Inorganic-ion transport and metabolism" =
    "Ion transport",

  "Function unknown" =
    "Function unknown"
)

biology_order <- c(
  "Carbohydrate metabolism",
  "Energy production",
  "Lipid metabolism",
  "Amino-acid metabolism",
  "Coenzyme metabolism",
  "Secondary metabolism",
  "Ion transport",
  "Translation",
  "Transcription",
  "RNA processing",
  "Chromatin",
  "Protein turnover",
  "Intracellular trafficking",
  "Cytoskeleton",
  "Cell wall and membrane",
  "Cell cycle"
)

# ============================================================
# Enrichment plotting data
# ============================================================

enrichment_plot_data <- enrichment %>%
  filter(
    padj < alpha,
    cog_name != "Function unknown"
  ) %>%
  mutate(
    label = recode(
      cog_name,
      !!!label_dictionary,
      .default = cog_name
    ),

    direction = factor(
      direction,
      levels = c("Induced", "Repressed")
    ),

    label = factor(
      label,
      levels = rev(biology_order)
    ),

    log2_odds_ratio = log2(odds_ratio)
  ) %>%
  filter(label %in% biology_order)

shared_enrichment_limit <- max(
  enrichment_plot_data$log2_odds_ratio,
  na.rm = TRUE
)

shared_enrichment_limit <- ceiling(
  shared_enrichment_limit * 10
) / 10

# ============================================================
# Function for species enrichment panels
# ============================================================

make_species_panel <- function(data, species_name, colour) {

  species_data <- data %>%
    filter(organism == species_name)

  ggplot(
    species_data,
    aes(
      x = log2_odds_ratio,
      y = label
    )
  ) +
    geom_col(
      fill = colour,
      width = 0.72
    ) +
    facet_grid(
      rows = vars(direction),
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
    scale_x_continuous(
      limits = c(0, shared_enrichment_limit),
      expand = expansion(
        mult = c(0, 0.04)
      )
    ) +
    labs(
      x = expression(log[2]~"enrichment odds ratio"),
      y = NULL
    ) +
    theme_classic(
      base_size = 12
    ) +
    theme(
      strip.background = element_blank(),

      strip.text.y.left = element_text(
        angle = 0,
        face = "bold",
        size = 11
      ),

      axis.text.y = element_text(
        size = 10
      ),

      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),

      panel.spacing.y = unit(
        5,
        "mm"
      ),

      plot.margin = margin(
        10,
        15,
        10,
        10
      )
    )
}

# ============================================================
# Panel A — Fusarium
# ============================================================

panel_a <- make_species_panel(
  enrichment_plot_data,
  "Fusarium",
  fusarium_colour
)

# ============================================================
# Panel B — Ophiostoma
# ============================================================

panel_b <- make_species_panel(
  enrichment_plot_data,
  "Ophiostoma",
  ophiostoma_colour
)

# ============================================================
# Panel C data — species × COG interactions
# ============================================================

interaction_plot_data <- interaction %>%
  filter(
    interaction_padj < alpha,
    cog_name != "Function unknown",
    fusarium_response_hits >= 20 |
      ophiostoma_response_hits >= 20
  ) %>%
  mutate(
    label = recode(
      cog_name,
      !!!label_dictionary,
      .default = cog_name
    ),

    response = factor(
      response,
      levels = c("Induced", "Repressed")
    ),

    log2_interaction = log2(
      interaction_odds_ratio
    ),

    log2_ci_lower = log2(
      interaction_ci_lower
    ),

    log2_ci_upper = log2(
      interaction_ci_upper
    ),

    supported_in = if_else(
      interaction_odds_ratio < 1,
      "Fusarium",
      "Ophiostoma"
    ),

    label = factor(
      label,
      levels = rev(biology_order)
    )
  ) %>%
  filter(label %in% biology_order)

interaction_limit <- max(
  abs(
    c(
      interaction_plot_data$log2_ci_lower,
      interaction_plot_data$log2_ci_upper
    )
  ),
  na.rm = TRUE
)

interaction_limit <- ceiling(
  interaction_limit * 2
) / 2

# ============================================================
# Panel C — interaction forest plot
# ============================================================

panel_c <- ggplot(
  interaction_plot_data,
  aes(
    x = log2_interaction,
    y = label
  )
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.45,
    linetype = "dashed",
    colour = "grey65"
  ) +
  geom_errorbarh(
    aes(
      xmin = log2_ci_lower,
      xmax = log2_ci_upper,
      colour = supported_in
    ),
    height = 0,
    linewidth = 0.7
  ) +
  geom_point(
    aes(
      colour = supported_in
    ),
    size = 3.2
  ) +
  facet_wrap(
    vars(response),
    ncol = 1,
    scales = "free_y",
    strip.position = "top"
  ) +
  scale_colour_manual(
    values = c(
      Fusarium = fusarium_colour,
      Ophiostoma = ophiostoma_colour
    ),
    breaks = c(
      "Fusarium",
      "Ophiostoma"
    ),
    labels = c(
      "Stronger in Fusarium",
      "Stronger in Ophiostoma"
    ),
    name = NULL
  ) +
  scale_x_continuous(
    limits = c(
      -interaction_limit,
      interaction_limit
    ),
    expand = expansion(
      mult = c(0.03, 0.03)
    )
  ) +
  labs(
    x = expression(
      log[2]~"species × COG interaction odds ratio"
    ),
    y = NULL
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    strip.background = element_blank(),

    strip.text.x = element_text(
      face = "bold",
      size = 11,
      hjust = 0
    ),

    axis.text.y = element_text(
      size = 10
    ),

    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),

    legend.position = "bottom",

    panel.spacing.y = unit(
      5,
      "mm"
    ),

    plot.margin = margin(
      10,
      15,
      10,
      10
    )
  )

# ============================================================
# Export separate panels
# ============================================================

ggsave(
  file.path(
    output_dir,
    "Figure4A_Fusarium_COG_enrichment.png"
  ),
  panel_a,
  width = 6.5,
  height = 7,
  dpi = 600
)

ggsave(
  file.path(
    output_dir,
    "Figure4A_Fusarium_COG_enrichment.pdf"
  ),
  panel_a,
  width = 6.5,
  height = 7
)

ggsave(
  file.path(
    output_dir,
    "Figure4B_Ophiostoma_COG_enrichment.png"
  ),
  panel_b,
  width = 6.5,
  height = 7,
  dpi = 600
)

ggsave(
  file.path(
    output_dir,
    "Figure4B_Ophiostoma_COG_enrichment.pdf"
  ),
  panel_b,
  width = 6.5,
  height = 7
)

ggsave(
  file.path(
    output_dir,
    "Figure4C_species_COG_interaction.png"
  ),
  panel_c,
  width = 7,
  height = 8,
  dpi = 600
)

ggsave(
  file.path(
    output_dir,
    "Figure4C_species_COG_interaction.pdf"
  ),
  panel_c,
  width = 7,
  height = 8
)

write_tsv(
  enrichment_plot_data,
  file.path(
    output_dir,
    "Figure4AB_enrichment_source_data.tsv"
  )
)

write_tsv(
  interaction_plot_data,
  file.path(
    output_dir,
    "Figure4C_interaction_source_data.tsv"
  )
)

cat(
  "\n",
  "============================================================\n",
  "FIGURE 4 PANELS COMPLETE\n",
  "============================================================\n\n",
  normalizePath(output_dir),
  "\n\n",
  sep = ""
)
