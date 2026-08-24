#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggdist)
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

point_sample_size <- 450L

point_data <- lfc |>
  group_by(organism) |>
  group_modify(
    ~ slice_sample(
      .x,
      n = min(point_sample_size, nrow(.x))
    )
  ) |>
  ungroup()

summary_data <- lfc |>
  group_by(organism) |>
  summarise(
    significant_genes = n(),
    median_abs_lfc = median(abs_log2_fold_change, na.rm = TRUE),
    q25_abs_lfc = quantile(
      abs_log2_fold_change, 0.25, na.rm = TRUE, names = FALSE
    ),
    q75_abs_lfc = quantile(
      abs_log2_fold_change, 0.75, na.rm = TRUE, names = FALSE
    ),
    p99_abs_lfc = quantile(
      abs_log2_fold_change, 0.99, na.rm = TRUE, names = FALSE
    ),
    maximum_abs_lfc = max(abs_log2_fold_change, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(
  summary_data,
  file.path(output_dir, "raincloud_distribution_summary.tsv")
)

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

p_final <- ggplot(
  lfc,
  aes(
    x = abs_log2_fold_change,
    y = organism_y,
    fill = organism
  )
) +

  # Cloud: clearly above the boxplot
  stat_halfeye(
    adjust = 0.7,
    width = 0.34,
    justification = -0.68,
    .width = 0,
    point_colour = NA,
    slab_alpha = 0.82,
    slab_colour = "#303030",
    slab_linewidth = 0.35
  ) +

  # Boxplot: centred on the species line
  geom_boxplot(
    aes(group = organism),
    width = 0.10,
    outlier.shape = NA,
    alpha = 0.95,
    linewidth = 0.38,
    colour = "#303030",
    position = position_nudge(y = -0.02)
  ) +

  # Median point
  stat_summary(
    aes(group = organism),
    fun = median,
    geom = "point",
    shape = 21,
    size = 2.15,
    stroke = 0.42,
    fill = "white",
    colour = "#202020",
    position = position_nudge(y = -0.02)
  ) +

  # Rain: shifted clearly below the boxplot
  geom_point(
    data = point_data,
    aes(
      x = abs_log2_fold_change,
      y = organism_y - 0.20,
      colour = organism
    ),
    inherit.aes = FALSE,
    position = position_jitter(
      width = 0,
      height = 0.034,
      seed = 114
    ),
    size = 0.72,
    alpha = 0.22,
    stroke = 0
  ) +

  scale_fill_manual(values = organism_fill) +
  scale_colour_manual(values = organism_fill) +

  scale_y_continuous(
    breaks = c(1, 2),
    labels = c(
      expression(italic("Fusarium")~"cf. salinense"),
      expression(italic("Ophiostoma")~"novo-ulmi")
    ),
    limits = c(0.60, 2.52),
    expand = expansion(mult = c(0, 0))
  ) +

  scale_x_continuous(
    breaks = pretty_breaks(n = 5),
    expand = expansion(mult = c(0.01, 0.04))
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
    axis.text.y = element_text(size = 10),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank(),
    plot.margin = margin(8, 10, 8, 8)
  )

ggsave(
  file.path(output_dir, "figure3D_raincloud_final.pdf"),
  p_final,
  width = 5.6,
  height = 3.15,
  units = "in",
  device = cairo_pdf
)

ggsave(
  file.path(output_dir, "figure3D_raincloud_final.png"),
  p_final,
  width = 5.6,
  height = 3.15,
  units = "in",
  dpi = 600,
  bg = "white"
)

message("\nGenerated final raincloud panel:")
message("  figure3D_raincloud_final.pdf")
message("  figure3D_raincloud_final.png")
message("\nFiles written to:")
message("  ", output_dir)
message("\nDisplayed x-axis maximum: ", shared_p99)
message("\nDistribution summary:")
print(summary_data, n = Inf)
