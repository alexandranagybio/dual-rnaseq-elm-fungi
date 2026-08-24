#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(ggplot2)
})

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
  "figure4_cog_enrichment"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

enrichment <- read_tsv(
  input_file,
  show_col_types = FALSE
)

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

plot_data <- enrichment %>%
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

    organism = factor(
      organism,
      levels = c("Fusarium", "Ophiostoma")
    ),

    direction = factor(
      direction,
      levels = c("Induced", "Repressed")
    ),

    log2_or = log2(odds_ratio)
  ) %>%
  filter(label %in% biology_order) %>%
  select(
    direction,
    organism,
    label,
    log2_or
  ) %>%
  complete(
    direction,
    organism,
    label,
    fill = list(log2_or = 0)
  ) %>%
  mutate(
    label = factor(
      label,
      levels = rev(biology_order)
    )
  )

cog_barplot <- ggplot(
  plot_data,
  aes(
    x = log2_or,
    y = label,
    fill = organism
  )
) +
  geom_col(
    position = position_dodge2(
      width = 0.82,
      preserve = "single"
    ),
    width = 0.72
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
    name = NULL
  ) +
  scale_x_continuous(
    limits = c(0, NA),
    expand = expansion(
      mult = c(0, 0.05)
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

    strip.text = element_text(
      face = "bold",
      size = 12
    ),

    axis.text.y = element_text(
      size = 10
    ),

    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),

    legend.position = "bottom",

    panel.spacing = unit(
      6,
      "mm"
    ),

    plot.margin = margin(
      10,
      15,
      10,
      10
    )
  )

ggsave(
  filename = file.path(
    output_dir,
    "Figure4A_COG_grouped_barplot.png"
  ),
  plot = cog_barplot,
  width = 8,
  height = 9,
  dpi = 600
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure4A_COG_grouped_barplot.pdf"
  ),
  plot = cog_barplot,
  width = 8,
  height = 9
)

cat(
  "\nCOG grouped bar plot complete\n\n",
  normalizePath(
    file.path(
      output_dir,
      "Figure4A_COG_grouped_barplot.png"
    )
  ),
  "\n\n",
  sep = ""
)
