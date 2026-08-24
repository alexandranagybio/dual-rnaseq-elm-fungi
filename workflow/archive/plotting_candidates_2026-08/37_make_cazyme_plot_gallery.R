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
# Input and output
# ==========================================================================

class_input <- file.path(
  "results",
  "publication",
  "cazyme_comparison",
  "plotting_tables",
  "cazyme_class_summary_long.tsv"
)

family_input <- file.path(
  "results",
  "publication",
  "cazyme_comparison",
  "plotting_tables",
  "cazyme_family_plotting_long.tsv"
)

output_dir <- file.path(
  "figures",
  "figure5_secretome_cazymes",
  "gallery"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ==========================================================================
# Colours
# ==========================================================================

fusarium_dark <- "#C86428"
fusarium_light <- "#E5A477"

ophiostoma_dark <- "#72508C"
ophiostoma_light <- "#A991BA"

neutral_dark <- "#343434"
neutral_light <- "#D9D9D9"

species_colours <- c(
  "Fusarium" = fusarium_dark,
  "Ophiostoma" = ophiostoma_dark
)

species_light_colours <- c(
  "Fusarium" = fusarium_light,
  "Ophiostoma" = ophiostoma_light
)

class_order <- c(
  "GH",
  "GT",
  "AA",
  "CBM",
  "CE",
  "PL"
)

species_order <- c(
  "Fusarium",
  "Ophiostoma"
)

# ==========================================================================
# Shared theme
# ==========================================================================

theme_publication <- function(base_size = 10) {
  theme_classic(
    base_size = base_size,
    base_family = "sans"
  ) +
    theme(
      axis.line = element_line(
        linewidth = 0.35,
        colour = "grey35"
      ),
      axis.ticks = element_line(
        linewidth = 0.3,
        colour = "grey40"
      ),
      axis.text = element_text(
        colour = "grey15"
      ),
      axis.title = element_text(
        colour = "grey15"
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        size = base_size,
        face = "bold",
        colour = "grey15"
      ),
      legend.title = element_text(
        size = base_size - 1
      ),
      legend.text = element_text(
        size = base_size - 1
      ),
      panel.spacing = unit(
        1.0,
        "lines"
      ),
      plot.margin = margin(
        8,
        8,
        8,
        8
      )
    )
}

# ==========================================================================
# Read and validate
# ==========================================================================

if (!file.exists(class_input)) {
  stop("Missing class table: ", class_input)
}

if (!file.exists(family_input)) {
  stop("Missing family table: ", family_input)
}

class_data <- read_tsv(
  class_input,
  show_col_types = FALSE,
  progress = FALSE
)

family_data <- read_tsv(
  family_input,
  show_col_types = FALSE,
  progress = FALSE
)

class_required <- c(
  "class",
  "species",
  "significant_genes",
  "up",
  "down",
  "secreted",
  "secreted_up",
  "secreted_down"
)

family_required <- c(
  "family",
  "class",
  "species",
  "significant_genes",
  "up",
  "down",
  "secreted",
  "secreted_up",
  "secreted_down",
  "mean_shrunk_lfc",
  "max_abs_shrunk_lfc"
)

missing_class <- setdiff(
  class_required,
  names(class_data)
)

missing_family <- setdiff(
  family_required,
  names(family_data)
)

if (length(missing_class) > 0) {
  stop(
    "Class table missing columns: ",
    paste(missing_class, collapse = ", ")
  )
}

if (length(missing_family) > 0) {
  stop(
    "Family table missing columns: ",
    paste(missing_family, collapse = ", ")
  )
}

# Add explicitly absent class combinations.
class_complete <- expand_grid(
  species = species_order,
  class = class_order
) %>%
  left_join(
    class_data,
    by = c(
      "species",
      "class"
    )
  ) %>%
  mutate(
    across(
      c(
        significant_genes,
        up,
        down,
        secreted,
        secreted_up,
        secreted_down
      ),
      ~ replace_na(.x, 0)
    ),
    species = factor(
      species,
      levels = species_order
    ),
    class = factor(
      class,
      levels = rev(class_order)
    )
  )

# ==========================================================================
# Candidate 1
# Class dumbbells: total versus secreted
# ==========================================================================

candidate_1_data <- class_complete %>%
  select(
    species,
    class,
    significant_genes,
    secreted
  )

candidate_1 <- ggplot(
  candidate_1_data,
  aes(
    y = class
  )
) +
  geom_segment(
    aes(
      x = secreted,
      xend = significant_genes,
      yend = class
    ),
    linewidth = 1.1,
    colour = "grey80",
    lineend = "round"
  ) +
  geom_point(
    aes(
      x = significant_genes
    ),
    size = 3.8,
    shape = 21,
    fill = "white",
    colour = neutral_dark,
    stroke = 0.8
  ) +
  geom_point(
    aes(
      x = secreted,
      colour = species
    ),
    size = 3.8
  ) +
  facet_wrap(
    vars(species),
    nrow = 1
  ) +
  scale_colour_manual(
    values = species_colours,
    guide = "none"
  ) +
  scale_x_continuous(
    expand = expansion(
      mult = c(0.03, 0.10)
    )
  ) +
  labs(
    x = "Genes",
    y = NULL
  ) +
  theme_publication() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# ==========================================================================
# Candidate 2
# Class Cleveland dots: induced versus repressed
# ==========================================================================

candidate_2_data <- class_complete %>%
  select(
    species,
    class,
    up,
    down
  ) %>%
  pivot_longer(
    cols = c(up, down),
    names_to = "direction",
    values_to = "genes"
  ) %>%
  mutate(
    direction = recode(
      direction,
      up = "Induced",
      down = "Repressed"
    ),
    direction = factor(
      direction,
      levels = c(
        "Induced",
        "Repressed"
      )
    )
  )

candidate_2 <- ggplot(
  candidate_2_data,
  aes(
    x = genes,
    y = class,
    shape = direction,
    colour = species
  )
) +
  geom_line(
    aes(
      group = interaction(
        species,
        class
      )
    ),
    colour = "grey84",
    linewidth = 0.8
  ) +
  geom_point(
    size = 3.8,
    stroke = 1
  ) +
  facet_wrap(
    vars(species),
    nrow = 1
  ) +
  scale_colour_manual(
    values = species_colours,
    guide = "none"
  ) +
  scale_shape_manual(
    values = c(
      "Induced" = 16,
      "Repressed" = 1
    )
  ) +
  labs(
    x = "Genes",
    y = NULL,
    shape = NULL
  ) +
  theme_publication() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "top"
  )

# ==========================================================================
# Select informative families
#
# Families score highly when they:
# - contain secreted genes;
# - show several significant genes;
# - differ in direction between species;
# - occur in only one species.
# ==========================================================================

family_scored <- family_data %>%
  mutate(
    secreted_fraction = if_else(
      significant_genes > 0,
      secreted / significant_genes,
      0
    ),
    directional_bias = abs(up - down),
    score =
      secreted * 4 +
      significant_genes * 1.5 +
      directional_bias * 2 +
      max_abs_shrunk_lfc
  ) %>%
  group_by(
    family,
    class
  ) %>%
  summarise(
    total_significant = sum(significant_genes),
    total_secreted = sum(secreted),
    maximum_score = max(score),
    species_present = n_distinct(species),
    .groups = "drop"
  ) %>%
  mutate(
    selection_score =
      maximum_score +
      if_else(
        species_present == 1,
        3,
        0
      )
  ) %>%
  arrange(
    desc(selection_score),
    desc(total_secreted),
    desc(total_significant),
    family
  ) %>%
  slice_head(n = 18)

selected_families <- family_scored$family

selected_family_data <- family_data %>%
  filter(
    family %in% selected_families
  ) %>%
  mutate(
    secreted_fraction = if_else(
      significant_genes > 0,
      secreted / significant_genes,
      0
    )
  )

family_order <- family_scored %>%
  arrange(selection_score) %>%
  pull(family)

# ==========================================================================
# Candidate 3
# Family bubble atlas
# ==========================================================================

candidate_3_data <- selected_family_data %>%
  select(
    family,
    class,
    species,
    up,
    down,
    secreted_up,
    secreted_down
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
    secreted = case_when(
      direction == "up" ~ secreted_up,
      direction == "down" ~ secreted_down,
      TRUE ~ 0
    ),
    secreted_fraction = if_else(
      genes > 0,
      secreted / genes,
      0
    ),
    column = paste(
      species,
      direction,
      sep = "_"
    ),
    column = factor(
      column,
      levels = c(
        "Fusarium_up",
        "Fusarium_down",
        "Ophiostoma_up",
        "Ophiostoma_down"
      ),
      labels = c(
        "Fusarium ↑",
        "Fusarium ↓",
        "Ophiostoma ↑",
        "Ophiostoma ↓"
      )
    ),
    family = factor(
      family,
      levels = family_order
    )
  ) %>%
  filter(
    genes > 0
  )

candidate_3 <- ggplot(
  candidate_3_data,
  aes(
    x = column,
    y = family
  )
) +
  geom_point(
    aes(
      size = genes,
      colour = species,
      alpha = secreted_fraction
    ),
    stroke = 0
  ) +
  scale_colour_manual(
    values = species_colours,
    guide = "none"
  ) +
  scale_size_area(
    max_size = 9,
    breaks = c(
      1,
      3,
      6
    ),
    name = "Genes"
  ) +
  scale_alpha_continuous(
    range = c(
      0.20,
      1
    ),
    limits = c(
      0,
      1
    ),
    labels = percent_format(
      accuracy = 1
    ),
    name = "Secreted"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_publication(
    base_size = 9
  ) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major = element_line(
      linewidth = 0.25,
      colour = "grey92"
    ),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    legend.position = "right"
  )

# ==========================================================================
# Candidate 4
# Family restrained heatmap
#
# Fill shows signed count:
# positive = induced
# negative = repressed
# Text shows the secreted count.
# ==========================================================================

candidate_4_data <- selected_family_data %>%
  mutate(
    net_direction = up - down,
    family = factor(
      family,
      levels = family_order
    ),
    species = factor(
      species,
      levels = species_order
    ),
    secreted_label = if_else(
      secreted > 0,
      as.character(secreted),
      ""
    )
  ) %>%
  complete(
    family,
    species,
    fill = list(
      net_direction = 0,
      secreted = 0,
      secreted_label = ""
    )
  )

heatmap_limit <- max(
  abs(candidate_4_data$net_direction),
  na.rm = TRUE
)

candidate_4 <- ggplot(
  candidate_4_data,
  aes(
    x = species,
    y = family,
    fill = net_direction
  )
) +
  geom_tile(
    width = 0.92,
    height = 0.92,
    colour = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(
      label = secreted_label
    ),
    size = 3,
    colour = "grey20"
  ) +
  scale_fill_gradient2(
    low = ophiostoma_dark,
    mid = "white",
    high = fusarium_dark,
    midpoint = 0,
    limits = c(
      -heatmap_limit,
      heatmap_limit
    ),
    name = "Up − down"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_publication(
    base_size = 9
  ) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

# ==========================================================================
# Candidate 5
# Secreted-family lollipop plot
#
# Focus only on secreted significant CAZymes.
# Positive = induced secreted genes.
# Negative = repressed secreted genes.
# ==========================================================================

candidate_5_data <- selected_family_data %>%
  select(
    family,
    class,
    species,
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
    signed_genes = if_else(
      direction == "secreted_down",
      -genes,
      genes
    ),
    family = factor(
      family,
      levels = family_order
    )
  ) %>%
  filter(
    genes > 0
  )

candidate_5 <- ggplot(
  candidate_5_data,
  aes(
    y = family,
    x = signed_genes,
    colour = species
  )
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.4,
    colour = "grey55"
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = signed_genes,
      yend = family
    ),
    linewidth = 0.8,
    alpha = 0.65
  ) +
  geom_point(
    size = 3.2
  ) +
  facet_wrap(
    vars(species),
    nrow = 1
  ) +
  scale_colour_manual(
    values = species_colours,
    guide = "none"
  ) +
  scale_x_continuous(
    labels = function(x) abs(x),
    expand = expansion(
      mult = c(
        0.08,
        0.08
      )
    )
  ) +
  labs(
    x = "Secreted genes",
    y = NULL
  ) +
  theme_publication(
    base_size = 9
  ) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# ==========================================================================
# Export individual candidates
# ==========================================================================

candidate_list <- list(
  "01_class_dumbbell" = candidate_1,
  "02_class_cleveland" = candidate_2,
  "03_family_bubble_atlas" = candidate_3,
  "04_family_heatmap" = candidate_4,
  "05_secreted_family_lollipop" = candidate_5
)

candidate_dimensions <- list(
  "01_class_dumbbell" = c(7.2, 4.0),
  "02_class_cleveland" = c(7.2, 4.0),
  "03_family_bubble_atlas" = c(7.4, 6.7),
  "04_family_heatmap" = c(5.7, 6.7),
  "05_secreted_family_lollipop" = c(7.4, 6.7)
)

for (candidate_name in names(candidate_list)) {
  candidate_plot <- candidate_list[[candidate_name]]
  dimensions <- candidate_dimensions[[candidate_name]]

  ggsave(
    filename = file.path(
      output_dir,
      paste0(
        candidate_name,
        ".pdf"
      )
    ),
    plot = candidate_plot,
    width = dimensions[1],
    height = dimensions[2],
    units = "in",
    device = cairo_pdf
  )

  ggsave(
    filename = file.path(
      output_dir,
      paste0(
        candidate_name,
        ".png"
      )
    ),
    plot = candidate_plot,
    width = dimensions[1],
    height = dimensions[2],
    units = "in",
    dpi = 450,
    bg = "white"
  )
}

# ==========================================================================
# Contact sheet
# ==========================================================================

contact_sheet <- (
  candidate_1 |
    candidate_2
) /
  (
    candidate_3 |
      candidate_4 |
      candidate_5
  ) +
  plot_layout(
    heights = c(
      0.8,
      1.45
    ),
    widths = c(
      1.35,
      1,
      1.35
    )
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

ggsave(
  filename = file.path(
    output_dir,
    "CAZyme_plot_gallery_contact_sheet.pdf"
  ),
  plot = contact_sheet,
  width = 15,
  height = 10,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    output_dir,
    "CAZyme_plot_gallery_contact_sheet.png"
  ),
  plot = contact_sheet,
  width = 15,
  height = 10,
  units = "in",
  dpi = 350,
  bg = "white"
)

# ==========================================================================
# Save selected families and summary
# ==========================================================================

write_tsv(
  family_scored,
  file.path(
    output_dir,
    "selected_family_ranking.tsv"
  )
)

cat("\n")
cat("============================================================\n")
cat("CAZYME PLOT GALLERY COMPLETE\n")
cat("============================================================\n\n")

cat(
  "Output directory:\n  ",
  output_dir,
  "\n\n",
  sep = ""
)

cat("Candidates:\n")

for (candidate_name in names(candidate_list)) {
  cat(
    "  ",
    candidate_name,
    ".png\n",
    sep = ""
  )
}

cat("\nContact sheet:\n")
cat(
  "  ",
  file.path(
    output_dir,
    "CAZyme_plot_gallery_contact_sheet.png"
  ),
  "\n\n",
  sep = ""
)

cat("Selected families:\n")

print(
  family_scored %>%
    select(
      family,
      class,
      total_significant,
      total_secreted,
      species_present,
      selection_score
    ),
  n = Inf,
  width = Inf
)
