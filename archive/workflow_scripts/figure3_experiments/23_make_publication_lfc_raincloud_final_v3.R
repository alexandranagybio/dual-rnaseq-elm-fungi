#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

set.seed(114)

repo_root <- normalizePath(".", mustWork = TRUE)

input_file <- file.path(
  repo_root,
  "results",
  "publication",
  "figure3",
  "figure3_lfc_significant.tsv"
)

output_dir <- file.path(
  repo_root,
  "results",
  "publication",
  "figure3",
  "raincloud_candidates"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop(
    "Missing input file:\n  ", input_file,
    "\nRun workflow/scripts/21_prepare_publication_lfc_distributions.R first."
  )
}

lfc <- read_tsv(
  input_file,
  show_col_types = FALSE,
  progress = FALSE
) |>
  filter(
    story == "Cross-species confrontation",
    significant,
    is.finite(abs_log2_fold_change)
  ) |>
  mutate(
    organism = factor(
      organism,
      levels = c("Fusarium", "Ophiostoma")
    ),
    organism_y = if_else(organism == "Fusarium", 1, 2)
  )

if (nrow(lfc) == 0) {
  stop("No significant cross-species confrontation genes found.")
}

# -------------------------------------------------------------------------
# Display settings
# -------------------------------------------------------------------------

point_sample_size <- 450L
density_adjust <- 1.0
density_height <- 0.24
density_offset <- 0.23
rain_offset <- 0.27

organism_fill <- c(
  "Fusarium" = "#DD8628",
  "Ophiostoma" = "#7A5795"
)

shared_p99 <- quantile(
  lfc$abs_log2_fold_change,
  probs = 0.99,
  na.rm = TRUE,
  names = FALSE
)

shared_p99 <- ceiling(shared_p99 * 4) / 4

# -------------------------------------------------------------------------
# Sample points for the decorative rain layer
# -------------------------------------------------------------------------

point_data <- lfc |>
  group_by(organism) |>
  group_modify(
    ~ slice_sample(
      .x,
      n = min(point_sample_size, nrow(.x))
    )
  ) |>
  ungroup() |>
  mutate(
    rain_y = organism_y - rain_offset
  )

# -------------------------------------------------------------------------
# Distribution summaries
#
# Box = Q1 to Q3
# Whiskers = 5th to 95th percentiles
# White point = median
# -------------------------------------------------------------------------

summary_data <- lfc |>
  group_by(organism, organism_y) |>
  summarise(
    significant_genes = n(),
    p05_abs_lfc = quantile(
      abs_log2_fold_change, 0.05, na.rm = TRUE, names = FALSE
    ),
    q25_abs_lfc = quantile(
      abs_log2_fold_change, 0.25, na.rm = TRUE, names = FALSE
    ),
    median_abs_lfc = median(
      abs_log2_fold_change, na.rm = TRUE
    ),
    q75_abs_lfc = quantile(
      abs_log2_fold_change, 0.75, na.rm = TRUE, names = FALSE
    ),
    p95_abs_lfc = quantile(
      abs_log2_fold_change, 0.95, na.rm = TRUE, names = FALSE
    ),
    p99_abs_lfc = quantile(
      abs_log2_fold_change, 0.99, na.rm = TRUE, names = FALSE
    ),
    maximum_abs_lfc = max(
      abs_log2_fold_change, na.rm = TRUE
    ),
    .groups = "drop"
  )

write_tsv(
  summary_data,
  file.path(output_dir, "raincloud_distribution_summary.tsv")
)

# -------------------------------------------------------------------------
# Manually construct smooth density polygons
#
# Very low-density tails are trimmed so the cloud outline does not continue
# as a nearly flat baseline across the complete panel.
# -------------------------------------------------------------------------

density_data <- lfc |>
  group_by(organism, organism_y) |>
  group_modify(function(.x, .y) {

    dens <- density(
      .x$abs_log2_fold_change,
      from = 0,
      to = max(shared_p99, max(.x$abs_log2_fold_change)),
      n = 1024,
      adjust = density_adjust,
      na.rm = TRUE
    )

    density_tbl <- tibble(
      x = dens$x,
      density = dens$y
    )

    density_cutoff <- max(density_tbl$density) * 0.0025

    density_tbl <- density_tbl |>
      filter(
        density >= density_cutoff,
        x <= shared_p99
      )

    if (nrow(density_tbl) < 2) {
      stop("Density trimming removed too many points for ", .y$organism)
    }

    baseline <- .y$organism_y + density_offset

    density_tbl |>
      mutate(
        density_scaled = density / max(density) * density_height,
        baseline = baseline,
        y = baseline + density_scaled
      )
  }) |>
  ungroup()

# -------------------------------------------------------------------------
# Final panel
# -------------------------------------------------------------------------

box_half_height <- 0.050
whisker_cap_half_height <- 0.035

p_final <- ggplot() +

  # Cloud
  geom_ribbon(
    data = density_data,
    aes(
      x = x,
      ymin = baseline,
      ymax = y,
      fill = organism,
      group = organism
    ),
    alpha = 0.82,
    colour = "#303030",
    linewidth = 0.32
  ) +

  # 5th–95th percentile whisker
  geom_segment(
    data = summary_data,
    aes(
      x = p05_abs_lfc,
      xend = p95_abs_lfc,
      y = organism_y,
      yend = organism_y
    ),
    linewidth = 0.42,
    colour = "#303030"
  ) +

  # Whisker caps
  geom_segment(
    data = summary_data,
    aes(
      x = p05_abs_lfc,
      xend = p05_abs_lfc,
      y = organism_y - whisker_cap_half_height,
      yend = organism_y + whisker_cap_half_height
    ),
    linewidth = 0.42,
    colour = "#303030"
  ) +

  geom_segment(
    data = summary_data,
    aes(
      x = p95_abs_lfc,
      xend = p95_abs_lfc,
      y = organism_y - whisker_cap_half_height,
      yend = organism_y + whisker_cap_half_height
    ),
    linewidth = 0.42,
    colour = "#303030"
  ) +

  # IQR box
  geom_rect(
    data = summary_data,
    aes(
      xmin = q25_abs_lfc,
      xmax = q75_abs_lfc,
      ymin = organism_y - box_half_height,
      ymax = organism_y + box_half_height,
      fill = organism
    ),
    alpha = 0.95,
    colour = "#303030",
    linewidth = 0.42
  ) +

  # Median line inside box
  geom_segment(
    data = summary_data,
    aes(
      x = median_abs_lfc,
      xend = median_abs_lfc,
      y = organism_y - box_half_height,
      yend = organism_y + box_half_height
    ),
    linewidth = 0.46,
    colour = "#303030"
  ) +

  # White median point
  geom_point(
    data = summary_data,
    aes(
      x = median_abs_lfc,
      y = organism_y
    ),
    shape = 21,
    size = 2.05,
    stroke = 0.40,
    fill = "white",
    colour = "#202020"
  ) +

  # Rain
  geom_point(
    data = point_data,
    aes(
      x = abs_log2_fold_change,
      y = rain_y,
      colour = organism
    ),
    position = position_jitter(
      width = 0,
      height = 0.030,
      seed = 114
    ),
    size = 0.55,
    alpha = 0.18,
    stroke = 0
  ) +

  scale_fill_manual(values = organism_fill) +
  scale_colour_manual(values = organism_fill) +

  scale_y_continuous(
    breaks = NULL,
    limits = c(0.56, 2.56),
    expand = expansion(mult = c(0, 0))
  ) +

  scale_x_continuous(
    breaks = pretty_breaks(n = 5),
    expand = expansion(add = c(0.035, 0.06))
  ) +

  coord_cartesian(
    xlim = c(0, shared_p99),
    clip = "off"
  ) +

  labs(
    x = expression("|log"[2] * " fold change|"),
    y = NULL
  ) +

  theme_classic(base_size = 10.5) +

  theme(
    axis.title.x = element_text(
      size = 10.5,
      margin = margin(t = 9)
    ),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "none",
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank(),
    plot.margin = margin(8, 10, 8, 8)
  )

ggsave(
  file.path(output_dir, "figure3D_raincloud_final_v3.pdf"),
  p_final,
  width = 5.0,
  height = 3.0,
  units = "in",
  device = cairo_pdf
)

ggsave(
  file.path(output_dir, "figure3D_raincloud_final_v3.png"),
  p_final,
  width = 5.0,
  height = 3.0,
  units = "in",
  dpi = 600,
  bg = "white"
)

message("\nGenerated refined raincloud panel:")
message("  figure3D_raincloud_final_v3.pdf")
message("  figure3D_raincloud_final_v3.png")
message("\nFiles written to:")
message("  ", output_dir)
message("\nDisplayed x-axis maximum: ", shared_p99)
message("\nBox: 25th–75th percentiles")
message("Whiskers: 5th–95th percentiles")
message("White point: median")
message("\nDistribution summary:")
print(summary_data, n = Inf)
