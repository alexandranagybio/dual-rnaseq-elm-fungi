#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

# ==========================================================================
# Figure 6b — functional signatures of Ophiostoma spatial response programs
#
# Visualization only.
# Uses existing significant enrichment results from script 39.
# ==========================================================================

ROOT <- normalizePath(".", mustWork = TRUE)

binary_file <- file.path(
  ROOT,
  "results/publication/ophiostoma_spatial_functional_enrichment",
  "spatial_secretome_cazyme_enrichment.tsv"
)

cog_file <- file.path(
  ROOT,
  "results/publication/ophiostoma_spatial_functional_enrichment",
  "spatial_cog_enrichment_significant.tsv"
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
  "Figure6b_spatial_functional_signature.pdf"
)

output_png <- file.path(
  output_dir,
  "Figure6b_spatial_functional_signature.png"
)

# ==========================================================================
# Display settings
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

cog_labels <- c(
  A = "RNA processing",
  B = "Chromatin",
  E = "Amino-acid metabolism",
  G = "Carbohydrate metabolism",
  I = "Lipid metabolism",
  J = "Translation",
  L = "Replication / repair",
  Q = "Secondary metabolism"
)

enriched_colour <- "#C56F78"
depleted_colour <- "#5D86A8"
neutral_colour <- "#666666"

# ==========================================================================
# Read
# ==========================================================================

for (f in c(binary_file, cog_file)) {
  if (!file.exists(f)) {
    stop("Missing input: ", f)
  }
}

binary <- read_tsv(
  binary_file,
  show_col_types = FALSE,
  progress = FALSE
)

cog <- read_tsv(
  cog_file,
  show_col_types = FALSE,
  progress = FALSE
)

# ==========================================================================
# Prepare significant functional labels
# ==========================================================================

binary_sig <- binary %>%
  filter(significant) %>%
  transmute(
    response_group,
    function_label = case_when(
      feature == "Secreted protein" ~ "Secreted proteins",
      feature == "High-confidence CAZyme" ~ "CAZymes",
      TRUE ~ feature
    ),
    direction = enrichment_direction,
    odds_ratio,
    padj
  )

cog_sig <- cog %>%
  filter(significant) %>%
  transmute(
    response_group,
    function_label = unname(
      cog_labels[cog]
    ),
    direction = enrichment_direction,
    odds_ratio,
    padj
  )

sig <- bind_rows(
  binary_sig,
  cog_sig
) %>%
  mutate(
    side = case_when(
      direction == "Enriched" ~ 1,
      direction == "Depleted" ~ -1,
      TRUE ~ 0
    )
  )

# ==========================================================================
# Build one display row per significant function
# ==========================================================================

# Rank labels within each response × direction group so they can be stacked.
plot_data <- sig %>%
  group_by(
    response_group,
    side
  ) %>%
  arrange(
    desc(abs(log2(odds_ratio))),
    function_label,
    .by_group = TRUE
  ) %>%
  mutate(
    item_index = row_number(),
    n_items = n()
  ) %>%
  ungroup()

# Base y position for each response program.
row_y <- tibble(
  response_group = response_order,
  base_y = rev(seq_along(response_order))
)

plot_data <- plot_data %>%
  left_join(
    row_y,
    by = "response_group"
  ) %>%
  mutate(
    # Small vertical offsets within each program.
    y = base_y +
      (
        item_index -
          (n_items + 1) / 2
      ) * 0.13,

    x = if_else(
      side < 0,
      -1,
      1
    ),

    hjust = if_else(
      side < 0,
      1,
      0
    )
  )

# ==========================================================================
# Reaction-zone note
# ==========================================================================

reaction_has_sig <- any(
  sig$response_group ==
    "Reaction-zone specific"
)

reaction_note <- tibble(
  response_group = "Reaction-zone specific",
  base_y = 4,
  x = 0,
  label = if (
    reaction_has_sig
  ) {
    ""
  } else {
    "No significant enrichment/depletion"
  }
)

# ==========================================================================
# Plot
# ==========================================================================

p <- ggplot() +

  # Central divider.
  geom_vline(
    xintercept = 0,
    linewidth = 0.40,
    colour = "grey78"
  ) +

  # Section headings.
  annotate(
    "text",
    x = -1,
    y = 4.72,
    label = "Depleted",
    family = "Nimbus Sans",
    fontface = "bold",
    size = 3.35,
    hjust = 1
  ) +

  annotate(
    "text",
    x = 1,
    y = 4.72,
    label = "Enriched",
    family = "Nimbus Sans",
    fontface = "bold",
    size = 3.35,
    hjust = 0
  ) +

  # Significant functional signatures.
  geom_point(
    data = plot_data,
    aes(
      x = x,
      y = y,
      colour = direction
    ),
    size = 2.4
  ) +

  geom_text(
    data = plot_data,
    aes(
      x = x +
        if_else(side < 0, -0.07, 0.07),
      y = y,
      label = function_label,
      colour = direction,
      hjust = hjust
    ),
    family = "Nimbus Sans",
    size = 3.1,
    show.legend = FALSE
  ) +

  # Explicit reaction-zone result.
  geom_text(
    data = reaction_note,
    aes(
      x = x,
      y = base_y,
      label = label
    ),
    family = "Nimbus Sans",
    size = 3.0,
    colour = neutral_colour
  ) +

  # Program labels.
  geom_text(
    data = row_y,
    aes(
      x = -2.30,
      y = base_y,
      label = response_labels[response_group]
    ),
    family = "Nimbus Sans",
    size = 3.1,
    hjust = 1,
    colour = "#222222"
  ) +

  scale_colour_manual(
    values = c(
      Enriched = enriched_colour,
      Depleted = depleted_colour
    ),
    guide = "none"
  ) +

  coord_cartesian(
    xlim = c(
      -2.45,
      2.45
    ),
    ylim = c(
      0.45,
      4.95
    ),
    clip = "off"
  ) +

  labs(
    x = NULL,
    y = NULL,
    tag = "b"
  ) +

  theme_void(
    base_family = "Nimbus Sans",
    base_size = 10
  ) +

  theme(
    plot.tag = element_text(
      face = "bold",
      size = 10,
      colour = "#222222"
    ),

    plot.tag.position = c(
      0.005,
      1.00
    ),

    plot.margin = margin(
      8,
      12,
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
  height = 92,
  units = "mm",
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  output_png,
  p,
  width = 180,
  height = 92,
  units = "mm",
  dpi = 600,
  bg = "white"
)

message("")
message("==============================================")
message("FIGURE 6B FUNCTIONAL SIGNATURE COMPLETE")
message("==============================================")
message("")
message("PDF: ", output_pdf)
message("PNG: ", output_png)
message("")
message("PASS: visualization only; no re-analysis.")
