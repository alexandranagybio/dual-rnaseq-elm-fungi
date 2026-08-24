#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ==========================================================================
# Figure 6b — extracellular-function enrichment across spatial programs
#
# Visualization only.
# Uses existing Fisher-test results from script 39.
# ==========================================================================

ROOT <- normalizePath(".", mustWork = TRUE)

input_file <- file.path(
  ROOT,
  "results/publication/ophiostoma_spatial_functional_enrichment",
  "spatial_secretome_cazyme_enrichment.tsv"
)

output_dir <- file.path(
  ROOT,
  "figures/publication/final"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_pdf <- file.path(
  output_dir,
  "Figure6b_spatial_extracellular_enrichment.pdf"
)

output_png <- file.path(
  output_dir,
  "Figure6b_spatial_extracellular_enrichment.png"
)

# ==========================================================================
# Display order / labels
# ==========================================================================

response_order <- c(
  "Reaction-zone specific",
  "Plate-wide confrontation response",
  "Complex spatial response",
  "Non-contact-region specific"
)

response_labels <- c(
  "Reaction-zone specific" =
    "Reaction-zone specific",

  "Plate-wide confrontation response" =
    "Plate-wide response",

  "Complex spatial response" =
    "Multi-region differential response",

  "Non-contact-region specific" =
    "Non-contact specific"
)

feature_labels <- c(
  "Secreted protein" = "Secreted proteins",
  "High-confidence CAZyme" = "CAZymes"
)

enriched_colour <- "#C56F78"
depleted_colour <- "#5D86A8"
neutral_colour <- "#888888"

# ==========================================================================
# Read / validate
# ==========================================================================

if (!file.exists(input_file)) {
  stop("Missing input: ", input_file)
}

dat <- read_tsv(
  input_file,
  show_col_types = FALSE,
  progress = FALSE
)

required_columns <- c(
  "response_group",
  "feature",
  "odds_ratio",
  "padj",
  "significant",
  "enrichment_direction"
)

missing_columns <- setdiff(
  required_columns,
  names(dat)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# ==========================================================================
# Prepare
# ==========================================================================

plot_data <- dat %>%
  filter(
    response_group %in% response_order,
    feature %in% names(feature_labels)
  ) %>%
  mutate(
    response_group = factor(
      response_group,
      levels = rev(response_order)
    ),

    response_label = factor(
      response_labels[as.character(response_group)],
      levels = rev(unname(response_labels[response_order]))
    ),

    feature_label = factor(
      feature_labels[feature],
      levels = c(
        "Secreted proteins",
        "CAZymes"
      )
    ),

    log2_or = log2(odds_ratio),

    state = case_when(
      significant & odds_ratio > 1 ~ "Enriched",
      significant & odds_ratio < 1 ~ "Depleted",
      TRUE ~ "Not significant"
    )
  )

x_limit <- max(
  abs(plot_data$log2_or),
  na.rm = TRUE
)

x_limit <- ceiling(
  x_limit * 4
) / 4

if (x_limit < 1) {
  x_limit <- 1
}

# ==========================================================================
# Plot
# ==========================================================================

p <- ggplot(
  plot_data,
  aes(
    x = log2_or,
    y = response_label
  )
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.40,
    linetype = "dashed",
    colour = "grey75"
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = log2_or,
      yend = response_label,
      colour = state
    ),
    linewidth = 0.70
  ) +
  geom_point(
    aes(
      colour = state,
      fill = state
    ),
    shape = 21,
    size = 3.8,
    stroke = 0.45
  ) +
  facet_wrap(
    ~ feature_label,
    nrow = 1
  ) +
  scale_colour_manual(
    values = c(
      "Enriched" = enriched_colour,
      "Depleted" = depleted_colour,
      "Not significant" = neutral_colour
    ),
    breaks = c(
      "Enriched",
      "Depleted",
      "Not significant"
    ),
    name = NULL
  ) +
  scale_fill_manual(
    values = c(
      "Enriched" = enriched_colour,
      "Depleted" = depleted_colour,
      "Not significant" = "white"
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(
      -x_limit,
      x_limit
    ),
    expand = expansion(
      mult = c(
        0.06,
        0.06
      )
    )
  ) +
  labs(
    x = expression(
      log[2]~"odds ratio"
    ),
    y = NULL,
    tag = "b"
  ) +
  theme_classic(
    base_size = 10,
    base_family = "Nimbus Sans"
  ) +
  theme(
    strip.background = element_blank(),

    strip.text = element_text(
      face = "bold",
      size = 9.5,
      margin = margin(
        b = 5
      )
    ),

    axis.text.x = element_text(
      size = 9
    ),

    axis.text.y = element_text(
      size = 9
    ),

    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),

    panel.spacing.x = unit(
      10,
      "mm"
    ),

    legend.position = "bottom",

    legend.text = element_text(
      size = 8.5
    ),

    plot.tag = element_text(
      face = "bold",
      size = 10
    ),

    plot.tag.position = c(
      0.005,
      1.02
    ),

    plot.margin = margin(
      8,
      10,
      6,
      8
    )
  )

# ==========================================================================
# Export
# ==========================================================================

ggsave(
  output_pdf,
  p,
  width = 180,
  height = 82,
  units = "mm",
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  output_png,
  p,
  width = 180,
  height = 82,
  units = "mm",
  dpi = 600,
  bg = "white"
)

message("")
message("==============================================")
message("FIGURE 6B EXTRACELLULAR ENRICHMENT COMPLETE")
message("==============================================")
message("")
message("PDF: ", output_pdf)
message("PNG: ", output_png)
message("")
message("PASS: visualization only; no re-analysis.")
