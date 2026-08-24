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

input_file <- file.path(
  "results",
  "publication",
  "cog_enrichment",
  "cog_enrichment_complete.tsv"
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
# Read enrichment results
# ============================================================

enrichment <- read_tsv(
  input_file,
  show_col_types = FALSE
)

# ============================================================
# Validate input
# ============================================================

required_columns <- c(
  "organism",
  "direction",
  "cog_name",
  "odds_ratio",
  "padj"
)

missing_columns <- setdiff(
  required_columns,
  names(enrichment)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# ============================================================
# Short category labels
# ============================================================

label_dictionary <- c(
  "Carbohydrate transport and metabolism" =
    "Carbohydrate",

  "Energy production and conversion" =
    "Energy",

  "Lipid transport and metabolism" =
    "Lipid",

  "Amino-acid transport and metabolism" =
    "Amino acid",

  "Coenzyme transport and metabolism" =
    "Coenzyme",

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
    "Trafficking",

  "Cytoskeleton" =
    "Cytoskeleton",

  "Cell wall, membrane and envelope biogenesis" =
    "Cell wall",

  "Cell-cycle control and chromosome partitioning" =
    "Cell cycle",

  "Secondary-metabolite biosynthesis and transport" =
    "Secondary metabolism",

  "Inorganic-ion transport and metabolism" =
    "Ion transport"
)

# ============================================================
# Biological category order
# ============================================================

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
  "Ion transport"
)

# ============================================================
# Prepare plotting data
#
# Only significant categories are retained.
# No artificial zero rows are added.
# ============================================================

plot_data <- enrichment %>%
  filter(
    padj < alpha,
    cog_name != "Function unknown",
    is.finite(odds_ratio),
    odds_ratio > 0
  ) %>%
  mutate(
    label = recode(
      cog_name,
      !!!label_dictionary,
      .default = cog_name
    ),

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

    log2_odds_ratio = log2(odds_ratio)
  ) %>%
  filter(
    label %in% biology_order
  ) %>%
  mutate(
    label = factor(
      label,
      levels = rev(biology_order)
    )
  )

if (nrow(plot_data) == 0) {
  stop(
    "No significantly enriched COG categories were available for plotting."
  )
}

# ============================================================
# Plot
# ============================================================

cog_barplot <- ggplot(
  plot_data,
  aes(
    x = log2_odds_ratio,
    y = label,
    fill = organism
  )
) +
  geom_col(
    position = position_dodge2(
      width = 0.82,
      preserve = "single"
    ),
    width = 0.70
  ) +
  facet_wrap(
    ~ direction,
    ncol = 1,
    scales = "free_y",
    strip.position = "top"
  ) +
  scale_fill_manual(
    values = c(
      Fusarium = fusarium_colour,
      Ophiostoma = ophiostoma_colour
    ),
    breaks = c(
      "Fusarium",
      "Ophiostoma"
    ),
    labels = c(
      expression(italic("Fusarium")~"cf. salinense"),
      expression(italic("Ophiostoma novo-ulmi"))
    ),
    name = NULL
  ) +
  scale_x_continuous(
    limits = c(0, NA),
    breaks = seq(0, 2, 0.5),
    expand = expansion(
      mult = c(0, 0.04)
    )
  ) +
  labs(
    x = expression(log[2]~"enrichment odds ratio"),
    y = NULL
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    strip.background = element_blank(),

    strip.text = element_text(
      face = "bold",
      size = 13,
      colour = "black",
      margin = margin(
        b = 8
      )
    ),

    panel.grid.major.x = element_line(
      colour = "grey88",
      linewidth = 0.4
    ),

    panel.grid.major.y = element_line(
      colour = "grey93",
      linewidth = 0.35
    ),

    panel.grid.minor = element_blank(),

    axis.text.y = element_text(
      size = 10.5,
      colour = "black"
    ),

    axis.text.x = element_text(
      size = 10,
      colour = "black"
    ),

    axis.title.x = element_text(
      size = 12,
      margin = margin(
        t = 8
      )
    ),

    axis.ticks = element_blank(),

    legend.position = "bottom",

    legend.text = element_text(
      size = 10.5
    ),

    legend.key.width = unit(
      8,
      "mm"
    ),

    panel.spacing.y = unit(
      7,
      "mm"
    ),

    plot.margin = margin(
      t = 10,
      r = 15,
      b = 10,
      l = 10
    )
  )

# ============================================================
# Export
# ============================================================

png_file <- file.path(
  output_dir,
  "Figure4AB_COG_enrichment_barplot.png"
)

pdf_file <- file.path(
  output_dir,
  "Figure4AB_COG_enrichment_barplot.pdf"
)

source_file <- file.path(
  output_dir,
  "Figure4AB_COG_enrichment_barplot_source_data.tsv"
)

ggsave(
  filename = png_file,
  plot = cog_barplot,
  width = 7.5,
  height = 7,
  dpi = 600
)

ggsave(
  filename = pdf_file,
  plot = cog_barplot,
  width = 7.5,
  height = 7
)

write_tsv(
  plot_data,
  source_file
)

cat(
  "\n",
  "============================================================\n",
  "COG ENRICHMENT BAR PLOT COMPLETE\n",
  "============================================================\n\n",
  "PNG:\n",
  normalizePath(png_file),
  "\n\nPDF:\n",
  normalizePath(pdf_file),
  "\n\nSource data:\n",
  normalizePath(source_file),
  "\n\n",
  sep = ""
)
