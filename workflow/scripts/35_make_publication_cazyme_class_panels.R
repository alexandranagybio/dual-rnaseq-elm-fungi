#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ==========================================================================
# Configuration
# ==========================================================================

input_file <- file.path(
  "results",
  "publication",
  "cazyme_comparison",
  "plotting_tables",
  "cazyme_class_summary_long.tsv"
)

output_dir <- file.path(
  "figures",
  "figure5_secretome_cazymes"
)

canonical_output_dir <- file.path(
  "results",
  "publication",
  "cazyme_comparison",
  "plotting_tables"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  canonical_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_pdf <- file.path(
  output_dir,
  "Figure5_CAZyme_class_response.pdf"
)

output_png <- file.path(
  output_dir,
  "Figure5_CAZyme_class_response.png"
)

# Figure-local copy for convenience.
output_table <- file.path(
  output_dir,
  "Figure5_CAZyme_class_response_plotting_table.tsv"
)

# Canonical table consumed by downstream inferential analyses.
canonical_output_table <- file.path(
  canonical_output_dir,
  "figure5_cazyme_class_response_plotting_table.tsv"
)

# Species colour language.
fusarium_dark <- "#C86428"
fusarium_light <- "#F2D8C5"

ophiostoma_dark <- "#72508C"
ophiostoma_light <- "#DDD0E5"

class_order <- c(
  "AA",
  "CBM",
  "CE",
  "GH",
  "GT",
  "PL"
)

species_order <- c(
  "Fusarium",
  "Ophiostoma"
)

# ==========================================================================
# Validation
# ==========================================================================

if (!file.exists(input_file)) {
  stop(
    "Missing input file: ",
    input_file
  )
}

required_columns <- c(
  "class",
  "species",
  "significant_genes",
  "up",
  "down",
  "secreted",
  "secreted_up",
  "secreted_down"
)

class_data <- read_tsv(
  input_file,
  show_col_types = FALSE,
  progress = FALSE
)

missing_columns <- setdiff(
  required_columns,
  names(class_data)
)

if (length(missing_columns) > 0) {
  stop(
    "Input table is missing columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

observed_species <- sort(
  unique(class_data$species)
)

expected_species <- sort(species_order)

if (!identical(
  observed_species,
  expected_species
)) {
  stop(
    "Unexpected species values. Observed: ",
    paste(
      observed_species,
      collapse = ", "
    )
  )
}

if (anyDuplicated(
  class_data[c("species", "class")]
) > 0) {
  stop(
    "Duplicated species × class rows detected."
  )
}

# Complete all expected species × class combinations so absent classes are
# represented explicitly as zero rather than being silently dropped.
class_data_complete <- tidyr::complete(
  class_data,
  species = species_order,
  class = class_order,
  fill = list(
    significant_genes = 0L,
    up = 0L,
    down = 0L,
    secreted = 0L,
    secreted_up = 0L,
    secreted_down = 0L
  )
)

validation <- class_data_complete %>%
  mutate(
    direction_sum = up + down,
    secreted_direction_sum = secreted_up + secreted_down
  )

if (any(
  validation$direction_sum !=
    validation$significant_genes
)) {
  stop(
    "One or more class rows fail: ",
    "up + down != significant_genes."
  )
}

if (any(
  validation$secreted_direction_sum !=
    validation$secreted
)) {
  stop(
    "One or more class rows fail: ",
    "secreted_up + secreted_down != secreted."
  )
}

if (any(
  validation$secreted >
    validation$significant_genes
)) {
  stop(
    "One or more class rows have more secreted ",
    "than significant genes."
  )
}

# ==========================================================================
# Long plotting table
# ==========================================================================

all_response <- class_data_complete %>%
  select(
    species,
    class,
    up,
    down
  ) %>%
  pivot_longer(
    cols = c(
      up,
      down
    ),
    names_to = "direction",
    values_to = "genes"
  ) %>%
  mutate(
    layer = "All significant",
    signed_genes = if_else(
      direction == "down",
      -genes,
      genes
    )
  )

secreted_response <- class_data_complete %>%
  select(
    species,
    class,
    secreted_up,
    secreted_down
  ) %>%
  pivot_longer(
    cols = c(
      secreted_up,
      secreted_down
    ),
    names_to = "direction",
    values_to = "genes"
  ) %>%
  mutate(
    direction = case_when(
      direction == "secreted_up" ~ "up",
      direction == "secreted_down" ~ "down",
      TRUE ~ NA_character_
    ),
    layer = "Secreted",
    signed_genes = if_else(
      direction == "down",
      -genes,
      genes
    )
  )

plot_data <- bind_rows(
  all_response,
  secreted_response
) %>%
  mutate(
    species = factor(
      species,
      levels = species_order
    ),
    class = factor(
      class,
      levels = rev(class_order)
    ),
    direction = factor(
      direction,
      levels = c(
        "down",
        "up"
      )
    ),
    layer = factor(
      layer,
      levels = c(
        "All significant",
        "Secreted"
      )
    )
  )

# Save both the figure-local copy and the canonical results-layer copy.
write_tsv(
  plot_data,
  output_table
)

write_tsv(
  plot_data,
  canonical_output_table
)

# ==========================================================================
# Shared scales
# ==========================================================================

maximum_value <- max(
  abs(plot_data$signed_genes),
  na.rm = TRUE
)

axis_limit <- ceiling(
  maximum_value / 10
) * 10

if (axis_limit == 0) {
  axis_limit <- 10
}

axis_breaks <- pretty(
  c(
    -axis_limit,
    axis_limit
  ),
  n = 5
)

# ==========================================================================
# Plot function
# ==========================================================================

make_species_panel <- function(
    species_name,
    panel_letter,
    light_colour,
    dark_colour) {

  species_data <- plot_data %>%
    filter(
      species == species_name
    )

  ggplot(
    species_data,
    aes(
      x = signed_genes,
      y = class
    )
  ) +
    geom_col(
      data = species_data %>%
        filter(
          layer == "All significant"
        ),
      aes(
        fill = "All significant"
      ),
      width = 0.70
    ) +
    geom_col(
      data = species_data %>%
        filter(
          layer == "Secreted"
        ),
      aes(
        fill = "Secreted"
      ),
      width = 0.38
    ) +
    geom_vline(
      xintercept = 0,
      linewidth = 0.35,
      colour = "#444444"
    ) +
    scale_x_continuous(
      limits = c(
        -axis_limit,
        axis_limit
      ),
      breaks = axis_breaks,
      labels = abs
    ) +
    scale_fill_manual(
      values = c(
        "All significant" = light_colour,
        "Secreted" = dark_colour
      ),
      breaks = c(
        "All significant",
        "Secreted"
      ),
      name = NULL
    ) +
    labs(
      x = "Genes",
      y = NULL,
      tag = panel_letter
    ) +
    theme_classic(
      base_size = 10
    ) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      plot.tag = element_text(
        face = "bold",
        size = 13
      ),
      axis.text.y = element_text(
        face = "bold"
      )
    )
}

fusarium_plot <- make_species_panel(
  "Fusarium",
  "A",
  fusarium_light,
  fusarium_dark
)

ophiostoma_plot <- make_species_panel(
  "Ophiostoma",
  "B",
  ophiostoma_light,
  ophiostoma_dark
)

combined_plot <- fusarium_plot /
  ophiostoma_plot +
  plot_layout(
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )

ggsave(
  output_pdf,
  combined_plot,
  width = 6.8,
  height = 7.2,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  output_png,
  combined_plot,
  width = 6.8,
  height = 7.2,
  units = "in",
  dpi = 600,
  bg = "white"
)

message("Wrote ", output_pdf)
message("Wrote ", output_png)
message("Wrote ", output_table)
message("Wrote canonical plot data: ", canonical_output_table)
