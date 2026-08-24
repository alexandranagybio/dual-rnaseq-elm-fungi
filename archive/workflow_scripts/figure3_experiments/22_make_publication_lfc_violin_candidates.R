#!/usr/bin/env Rscript

# 22_make_publication_lfc_violin_candidates.R
#
# Purpose:
#   Generate exploratory Figure 3D violin-plot candidates.
#
# The script creates four versions:
#   D1: signed log2FC, cross-species confrontation only
#   D2: absolute log2FC, cross-species confrontation only
#   D3: signed log2FC, all contrasts grouped by biological story
#   D4: signed log2FC with separate up/down half-violins approximated by facets
#
# Run after script 21:
#   Rscript workflow/scripts/22_make_publication_lfc_violin_candidates.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(stringr)
})

repo_root <- normalizePath(".", mustWork = TRUE)
input_dir <- file.path(repo_root, "results", "publication", "figure3")
output_dir <- file.path(input_dir, "violin_candidates")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(input_dir, "figure3_lfc_significant.tsv")
summary_file <- file.path(input_dir, "figure3_lfc_distribution_summary.tsv")

if (!file.exists(input_file)) {
  stop(
    "Missing input: ", input_file,
    "\nRun workflow/scripts/21_prepare_publication_lfc_distributions.R first."
  )
}

lfc <- read_tsv(
  input_file,
  show_col_types = FALSE,
  progress = FALSE
) |>
  mutate(
    contrast_display = str_replace_all(
      contrast_display,
      fixed(" | "),
      "\\n"
    )
  )

summary_tbl <- read_tsv(
  summary_file,
  show_col_types = FALSE,
  progress = FALSE
)

required_columns <- c(
  "organism", "contrast_id", "contrast_label", "story",
  "contrast_display", "gene_id", "log2_fold_change",
  "abs_log2_fold_change", "direction"
)
missing_columns <- setdiff(required_columns, names(lfc))
if (length(missing_columns) > 0) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

# -------------------------------------------------------------------------
# Visual constants
# Replace these with the exact manuscript palette when finalized.
# -------------------------------------------------------------------------

organism_fill <- c(
  "Fusarium" = "#D7832F",
  "Ophiostoma" = "#76528E"
)

direction_fill <- c(
  "Downregulated" = "#D7D0C5",
  "Upregulated" = "#5A5A5A"
)

base_theme <- theme_classic(base_size = 10) +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    axis.text.x = element_text(lineheight = 0.95),
    legend.title = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9),
    plot.margin = margin(7, 9, 7, 7)
  )

save_plot <- function(plot, stem, width, height) {
  ggsave(
    filename = file.path(output_dir, paste0(stem, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf
  )
  ggsave(
    filename = file.path(output_dir, paste0(stem, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )
}

add_n_labels <- function(data, y_position, value_prefix = "n = ") {
  data |>
    count(organism, contrast_display, name = "n") |>
    mutate(
      y = y_position,
      label = paste0(value_prefix, comma(n))
    )
}

# -------------------------------------------------------------------------
# D1 — Signed LFC, like-for-like cross-species comparison
#
# Biological question:
#   Among significant genes, how large and directionally distributed were
#   the responses to confrontation?
# -------------------------------------------------------------------------

cross_species <- lfc |>
  filter(story == "Cross-species confrontation") |>
  mutate(
    contrast_display = factor(
      contrast_display,
      levels = c(
        "Fusarium\nInteraction vs control",
        "Ophiostoma\nInteraction vs control"
      )
    )
  )

signed_limit <- max(abs(cross_species$log2_fold_change), na.rm = TRUE)
signed_limit <- ceiling(signed_limit * 2) / 2
signed_limit <- max(signed_limit, 1)

d1_labels <- add_n_labels(
  cross_species,
  y_position = signed_limit * 0.97
)

p_d1 <- ggplot(
  cross_species,
  aes(
    x = contrast_display,
    y = log2_fold_change,
    fill = organism
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    width = 0.78,
    linewidth = 0.35
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.35
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.35,
    linetype = "dashed"
  ) +
  geom_text(
    data = d1_labels,
    aes(x = contrast_display, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.7,
    vjust = 1
  ) +
  scale_fill_manual(values = organism_fill) +
  scale_y_continuous(
    limits = c(-signed_limit, signed_limit),
    breaks = pretty_breaks(n = 5),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = NULL,
    y = expression(log[2] * " fold change"),
    title = "Magnitude and direction of significant transcriptional changes",
    subtitle = "Violin width shows distribution shape; internal box shows median and interquartile range"
  ) +
  base_theme

save_plot(
  p_d1,
  "figure3D_candidate_1_signed_lfc_cross_species",
  width = 4.6,
  height = 3.8
)

# -------------------------------------------------------------------------
# D2 — Absolute LFC, like-for-like cross-species comparison
#
# Biological question:
#   Among genes that responded significantly, which fungus showed larger
#   expression changes irrespective of direction?
# -------------------------------------------------------------------------

abs_limit <- max(cross_species$abs_log2_fold_change, na.rm = TRUE)
abs_limit <- ceiling(abs_limit * 2) / 2

d2_labels <- add_n_labels(
  cross_species,
  y_position = abs_limit * 0.97
)

median_abs <- cross_species |>
  group_by(organism, contrast_display) |>
  summarise(
    median_abs_lfc = median(abs_log2_fold_change, na.rm = TRUE),
    .groups = "drop"
  )

p_d2 <- ggplot(
  cross_species,
  aes(
    x = contrast_display,
    y = abs_log2_fold_change,
    fill = organism
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    width = 0.78,
    linewidth = 0.35
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.35
  ) +
  geom_text(
    data = d2_labels,
    aes(x = contrast_display, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.7,
    vjust = 1
  ) +
  scale_fill_manual(values = organism_fill) +
  scale_y_continuous(
    limits = c(0, abs_limit),
    breaks = pretty_breaks(n = 5),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = expression("|log"[2] * " fold change|"),
    title = "Magnitude of significant transcriptional changes",
    subtitle = "Direction removed to compare response intensity"
  ) +
  base_theme

save_plot(
  p_d2,
  "figure3D_candidate_2_absolute_lfc_cross_species",
  width = 4.6,
  height = 3.8
)

# -------------------------------------------------------------------------
# D3 — Signed LFC across both biological stories
#
# This is useful for exploration, but it may be too information-dense for
# the final main figure.
# -------------------------------------------------------------------------

all_levels <- c(
  "Fusarium\nInteraction vs control",
  "Ophiostoma\nInteraction vs control",
  "Ophiostoma\nInteraction vs outside",
  "Ophiostoma\nOutside vs control"
)

all_contrasts <- lfc |>
  mutate(
    contrast_display = factor(contrast_display, levels = all_levels),
    story = factor(
      story,
      levels = c(
        "Cross-species confrontation",
        "Ophiostoma spatial response"
      )
    )
  )

all_signed_limit <- max(abs(all_contrasts$log2_fold_change), na.rm = TRUE)
all_signed_limit <- ceiling(all_signed_limit * 2) / 2

d3_labels <- all_contrasts |>
  count(story, organism, contrast_display, name = "n") |>
  mutate(
    y = all_signed_limit * 0.97,
    label = paste0("n = ", comma(n))
  )

p_d3 <- ggplot(
  all_contrasts,
  aes(
    x = contrast_display,
    y = log2_fold_change,
    fill = organism
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    width = 0.78,
    linewidth = 0.3
  ) +
  geom_boxplot(
    width = 0.11,
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.3
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.3,
    linetype = "dashed"
  ) +
  geom_text(
    data = d3_labels,
    aes(x = contrast_display, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.4,
    vjust = 1
  ) +
  facet_grid(
    cols = vars(story),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(values = organism_fill) +
  scale_y_continuous(
    limits = c(-all_signed_limit, all_signed_limit),
    breaks = pretty_breaks(n = 5),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = NULL,
    y = expression(log[2] * " fold change"),
    title = "Distributions of significant transcriptional changes"
  ) +
  base_theme +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    panel.spacing.x = unit(1.2, "lines")
  )

save_plot(
  p_d3,
  "figure3D_candidate_3_signed_lfc_all_contrasts",
  width = 7.2,
  height = 3.8
)

# -------------------------------------------------------------------------
# D4 — Direction-separated distributions
#
# Biological question:
#   Are up- and downregulated genes characterized by similar effect sizes?
#
# Absolute values are used after separating direction, which avoids plotting
# the same information on opposite sides of zero.
# -------------------------------------------------------------------------

direction_levels <- c("Downregulated", "Upregulated")

direction_data <- cross_species |>
  filter(direction %in% direction_levels) |>
  mutate(
    direction = factor(direction, levels = direction_levels),
    contrast_display = factor(
      contrast_display,
      levels = c(
        "Fusarium\nInteraction vs control",
        "Ophiostoma\nInteraction vs control"
      )
    )
  )

direction_labels <- direction_data |>
  count(direction, organism, contrast_display, name = "n") |>
  group_by(direction) |>
  mutate(
    y = max(direction_data$abs_log2_fold_change, na.rm = TRUE) * 0.97,
    label = paste0("n = ", comma(n))
  ) |>
  ungroup()

p_d4 <- ggplot(
  direction_data,
  aes(
    x = contrast_display,
    y = abs_log2_fold_change,
    fill = organism
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    width = 0.78,
    linewidth = 0.3
  ) +
  geom_boxplot(
    width = 0.11,
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.3
  ) +
  geom_text(
    data = direction_labels,
    aes(x = contrast_display, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.3,
    vjust = 1
  ) +
  facet_wrap(
    vars(direction),
    nrow = 1
  ) +
  scale_fill_manual(values = organism_fill) +
  scale_y_continuous(
    limits = c(0, max(direction_data$abs_log2_fold_change, na.rm = TRUE)),
    breaks = pretty_breaks(n = 5),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = expression("|log"[2] * " fold change|"),
    title = "Magnitude of up- and downregulated responses"
  ) +
  base_theme +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9)
  )

save_plot(
  p_d4,
  "figure3D_candidate_4_direction_separated",
  width = 7.0,
  height = 3.8
)

# -------------------------------------------------------------------------
# Compact numerical summary for interpretation
# -------------------------------------------------------------------------

interpretation_summary <- cross_species |>
  group_by(organism, direction) |>
  summarise(
    genes = n(),
    median_lfc = median(log2_fold_change, na.rm = TRUE),
    median_abs_lfc = median(abs_log2_fold_change, na.rm = TRUE),
    q25_abs_lfc = quantile(
      abs_log2_fold_change, 0.25, na.rm = TRUE, names = FALSE
    ),
    q75_abs_lfc = quantile(
      abs_log2_fold_change, 0.75, na.rm = TRUE, names = FALSE
    ),
    percent_abs_lfc_gt_1 = 100 * mean(abs_log2_fold_change > 1, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(
  interpretation_summary,
  file.path(output_dir, "figure3D_cross_species_distribution_summary.tsv")
)

message("Generated violin candidates:")
message("  D1 signed cross-species")
message("  D2 absolute cross-species")
message("  D3 signed all contrasts")
message("  D4 direction-separated cross-species")
message("\nFiles written to:")
message("  ", output_dir)
message("\nCross-species distribution summary:")
print(interpretation_summary, n = Inf)
