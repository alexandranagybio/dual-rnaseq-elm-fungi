#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
})

# ==========================================================================
# Figure 4
# Comparative COG enrichment heatmap
#
# Panel A
#   Comparative enrichment heatmap
#
# Panel B
#   Species × COG interaction plot
#
# Alexandra Nagy
# ==========================================================================

# ==========================================================================
# Configuration
# ==========================================================================

alpha_threshold <- 0.05

fusarium_colour <- "#D9792B"
ophiostoma_colour <- "#76519D"

neutral_colour <- "grey95"

heatmap_cell_size <- 0.95

# ==========================================================================
# Input files
# ==========================================================================

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

# ==========================================================================
# Output directories
# ==========================================================================

figure_dir <- file.path(
  "figures",
  "figure4_cog_enrichment"
)

table_dir <- file.path(
  "results",
  "publication",
  "figure4_cog_final"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ==========================================================================
# Load data
# ==========================================================================

message("Loading enrichment tables...")

enrichment <- read_tsv(
  enrichment_file,
  show_col_types = FALSE
)

interaction <- read_tsv(
  interaction_file,
  show_col_types = FALSE
)

# ==========================================================================
# Basic validation
# ==========================================================================

required_columns <- c(
  "organism",
  "direction",
  "cog",
  "cog_name",
  "odds_ratio",
  "padj"
)

missing_columns <- setdiff(
  required_columns,
  names(enrichment)
)

if(length(missing_columns) > 0){

  stop(
    "Missing columns:\n",
    paste(missing_columns, collapse="\n")
  )

}

# ==========================================================================
# Prepare enrichment table
# ==========================================================================

heatmap <- enrichment %>%

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
    )

  ) %>%

  mutate(

    column = case_when(

      organism == "Fusarium" &
        direction == "Induced"
      ~ "F ↑",

      organism == "Fusarium" &
        direction == "Repressed"
      ~ "F ↓",

      organism == "Ophiostoma" &
        direction == "Induced"
      ~ "O ↑",

      TRUE
      ~ "O ↓"

    )

  )

# ==========================================================================
# Short labels
# ==========================================================================

label_dictionary <- c(

  "Carbohydrate transport and metabolism" =
    "Carbohydrate",

  "Energy production and conversion" =
    "Energy",

  "Translation and ribosome biogenesis" =
    "Translation",

  "Amino-acid transport and metabolism" =
    "Amino acid",

  "Coenzyme transport and metabolism" =
    "Coenzyme",

  "Lipid transport and metabolism" =
    "Lipid",

  "RNA processing and modification" =
    "RNA processing",

  "Transcription" =
    "Transcription",

  "Chromatin structure and dynamics" =
    "Chromatin",

  "Post-translational modification and protein turnover" =
    "Protein turnover",

  "Intracellular trafficking and secretion" =
    "Trafficking",

  "Cell wall, membrane and envelope biogenesis" =
    "Cell wall",

  "Cell-cycle control and chromosome partitioning" =
    "Cell cycle",

  "Cytoskeleton" =
    "Cytoskeleton",

  "Secondary-metabolite biosynthesis and transport" =
    "Secondary metabolism",

  "Inorganic-ion transport and metabolism" =
    "Ion transport",

  "Function unknown" =
    "Unknown"

)

heatmap <- heatmap %>%

  mutate(

    label = recode(
      cog_name,
      !!!label_dictionary,
      .default = cog_name
    )

  )

# ==========================================================================
# Biological ordering
# ==========================================================================

biology_order <- c(

  "Carbohydrate",
  "Energy",
  "Lipid",
  "Amino acid",
  "Coenzyme",

  "Translation",
  "Transcription",
  "RNA processing",
  "Chromatin",

  "Protein turnover",
  "Trafficking",
  "Cytoskeleton",
  "Cell wall",
  "Cell cycle",

  "Secondary metabolism",
  "Ion transport",
  "Unknown"

)

heatmap <- heatmap %>%

  mutate(

    label = factor(
      label,
      levels = rev(biology_order)
    )

  )

# ==========================================================================
# Heatmap values
# ==========================================================================

heatmap <- heatmap %>%

  mutate(

    log2_or = log2(odds_ratio),

    significant = padj < alpha_threshold,

    fill_value = ifelse(
      significant,
      log2_or,
      0
    )

  )

# ==========================================================================
# Complete grid
# ==========================================================================

heatmap <- heatmap %>%

  complete(

    label,

    column = c(
      "F ↑",
      "F ↓",
      "O ↑",
      "O ↓"
    ),

    fill = list(

      fill_value = 0,

      significant = FALSE

    )

  )

# ==========================================================================
# Determine colour scale limits
# ==========================================================================

heatmap_limit <-

  max(
    abs(heatmap$fill_value),
    na.rm = TRUE
  )

message(
  "Heatmap range: ±",
  round(
    heatmap_limit,
    2
  )
)

# ==========================================================================
# Panel A
# ==========================================================================

panel_a <-

  ggplot(
    heatmap,
    aes(
      x = column,
      y = label,
      fill = fill_value
    )
  ) +

  geom_tile(
    colour = "white",
    linewidth = 0.5,
    width = 0.95,
    height = 0.95
  ) +

  geom_point(
    data = subset(
      heatmap,
      significant
    ),
    shape = 16,
    size = 1.3,
    colour = "black"
  ) +

  scale_fill_gradient2(

    low = fusarium_colour,
    mid = "white",
    high = ophiostoma_colour,

    midpoint = 0,

    limits = c(
      -heatmap_limit,
      heatmap_limit
    ),

    name = expression(log[2]("Odds ratio"))

  ) +

  labs(
    x = NULL,
    y = NULL
  ) +

  theme_classic(base_size = 11) +

  theme(

    axis.text.x = element_text(
      face = "bold",
      size = 11
    ),

    axis.text.y = element_text(
      size = 10
    ),

    axis.ticks = element_blank(),

    panel.border = element_blank(),

    legend.position = "right",

    plot.margin = margin(
      5,
      5,
      5,
      5
    )

  )

# ==========================================================================
# Panel B
# ==========================================================================

interaction_plot <-

  interaction %>%

  filter(

    interaction_padj < alpha_threshold,

    fusarium_response_hits >= 20 |
      ophiostoma_response_hits >= 20

  ) %>%

  mutate(

    label = recode(
      cog_name,
      !!!label_dictionary,
      .default = cog_name
    ),

    label = factor(
      label,
      levels = rev(biology_order)
    ),

    log2_interaction =
      log2(interaction_odds_ratio),

    colour = ifelse(

      interaction_odds_ratio > 1,

      "Ophiostoma",

      "Fusarium"

    )

  )

panel_b <-

  ggplot(

    interaction_plot,

    aes(

      x = log2_interaction,

      y = label

    )

  ) +

  geom_vline(

    xintercept = 0,

    colour = "grey70",

    linewidth = 0.4,

    linetype = 2

  ) +

  geom_errorbarh(

    aes(

      xmin = log2(interaction_ci_lower),

      xmax = log2(interaction_ci_upper)

    ),

    height = 0,

    linewidth = 0.5,

    colour = "grey55"

  ) +

  geom_point(

    aes(

      size = pmax(

        fusarium_response_hits,

        ophiostoma_response_hits

      ),

      colour = colour

    )

  ) +

  scale_colour_manual(

    values = c(

      Fusarium = fusarium_colour,

      Ophiostoma = ophiostoma_colour

    ),

    guide = "none"

  ) +

  scale_size_continuous(

    range = c(

      2,

      5

    ),

    name = "Genes"

  ) +

  labs(

    x = expression(log[2]("Interaction odds ratio")),

    y = NULL

  ) +

  theme_classic(base_size = 11) +

  theme(

    axis.text.y = element_blank(),

    axis.ticks.y = element_blank(),

    panel.border = element_blank(),

    plot.margin = margin(

      5,

      5,

      5,

      0

    )

  )

# ==========================================================================
# Assemble
# ==========================================================================

figure4 <-

  panel_a +

  panel_b +

  plot_layout(

    widths = c(

      1.15,

      1

    )

  )

# ==========================================================================
# Export
# ==========================================================================

ggsave(

  filename = file.path(

    figure_dir,

    "Figure4_COG_heatmap.pdf"

  ),

  plot = figure4,

  width = 11,

  height = 7,

  dpi = 600

)

ggsave(

  filename = file.path(

    figure_dir,

    "Figure4_COG_heatmap.png"

  ),

  plot = figure4,

  width = 11,

  height = 7,

  dpi = 600

)

write_tsv(

  heatmap,

  file.path(

    table_dir,

    "Figure4_heatmap_source_data.tsv"

  )

)

write_tsv(

  interaction_plot,

  file.path(

    table_dir,

    "Figure4_interaction_source_data.tsv"

  )

)

cat(

  "\n",

  "============================================================\n",

  "FIGURE 4 HEATMAP COMPLETE\n",

  "============================================================\n\n",

  "Figure written to:\n",

  normalizePath(figure_dir),

  "\n\n",

  sep = ""

)
