#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggdist)
  library(stringr)
  library(scales)
})

set.seed(114)

repo_root <- normalizePath(".", mustWork = TRUE)
input_file <- file.path(
  repo_root, "results", "publication", "figure3",
  "figure3_lfc_significant.tsv"
)
output_dir <- file.path(
  repo_root, "results", "publication", "figure3",
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
    )
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

raincloud_theme <- theme_classic(base_size = 10.5) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10.5, margin = margin(t = 8)),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 9),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 11.5, face = "bold"),
    plot.subtitle = element_text(size = 9.2, margin = margin(b = 8)),
    plot.caption = element_text(
      size = 7.7,
      hjust = 0,
      margin = margin(t = 8)
    ),
    plot.margin = margin(10, 14, 10, 10)
  )

make_raincloud <- function(data, points, x_limits = NULL,
                           x_trans = "identity",
                           title,
                           subtitle,
                           caption = NULL) {

  p <- ggplot(
    data,
    aes(
      x = abs_log2_fold_change,
      y = organism,
      fill = organism
    )
  ) +
    stat_halfeye(
      adjust = 0.7,
      width = 0.72,
      justification = -0.25,
      .width = 0,
      point_colour = NA,
      slab_alpha = 0.82,
      slab_colour = "#303030",
      slab_linewidth = 0.35
    ) +
    geom_boxplot(
      width = 0.13,
      outlier.shape = NA,
      alpha = 0.95,
      linewidth = 0.38,
      colour = "#303030",
      position = position_nudge(y = 0.01)
    ) +
    geom_point(
      data = points,
      aes(
        x = abs_log2_fold_change,
        y = organism,
        colour = organism
      ),
      inherit.aes = FALSE,
      position = position_jitter(
        width = 0,
        height = 0.075,
        seed = 114
      ),
      size = 0.75,
      alpha = 0.22,
      stroke = 0
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 21,
      size = 2.3,
      stroke = 0.45,
      fill = "white",
      colour = "#202020",
      position = position_nudge(y = 0.01)
    ) +
    scale_fill_manual(values = organism_fill) +
    scale_colour_manual(values = organism_fill) +
    scale_y_discrete(
      labels = c(
        "Fusarium" = expression(italic("Fusarium")~"cf. salinense"),
        "Ophiostoma" = expression(italic("Ophiostoma")~"novo-ulmi")
      )
    ) +
    labs(
      x = expression("|log"[2] * " fold change|"),
      title = title,
      subtitle = subtitle,
      caption = caption
    ) +
    raincloud_theme

  if (identical(x_trans, "identity")) {
    p <- p +
      scale_x_continuous(
        limits = x_limits,
        breaks = pretty_breaks(n = 5),
        expand = expansion(mult = c(0.01, 0.04))
      )
  } else if (identical(x_trans, "pseudo_log")) {
    p <- p +
      scale_x_continuous(
        trans = pseudo_log_trans(sigma = 0.2, base = 2),
        breaks = c(0, 0.25, 0.5, 1, 2, 4, 8),
        labels = label_number(accuracy = 0.01),
        expand = expansion(mult = c(0.01, 0.04))
      )
  }

  p
}

save_candidate <- function(plot, stem, width = 6.0, height = 3.5) {
  ggsave(
    file.path(output_dir, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf
  )

  ggsave(
    file.path(output_dir, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )
}

p_full <- make_raincloud(
  data = lfc,
  points = point_data,
  title = "Magnitude of significant transcriptional responses",
  subtitle = "Density and summary statistics use all significant genes",
  caption = paste0(
    "Points show a deterministic sample of up to ",
    point_sample_size,
    " genes per fungus; clouds and box summaries use the complete distributions."
  )
)

save_candidate(
  p_full,
  "figure3D_raincloud_candidate_1_full_range"
)

shared_p99 <- quantile(
  lfc$abs_log2_fold_change,
  probs = 0.99,
  na.rm = TRUE,
  names = FALSE
)
shared_p99 <- ceiling(shared_p99 * 4) / 4

outside_counts <- lfc |>
  group_by(organism) |>
  summarise(
    beyond_axis = sum(abs_log2_fold_change > shared_p99),
    .groups = "drop"
  )

outside_note <- paste(
  paste0(outside_counts$organism, ": ", outside_counts$beyond_axis),
  collapse = "; "
)

p_zoom <- make_raincloud(
  data = lfc,
  points = point_data,
  x_limits = c(0, shared_p99),
  title = "Magnitude of significant transcriptional responses",
  subtitle = "Ophiostoma responses were broader in effect magnitude",
  caption = paste0(
    "Linear axis displayed to the shared 99th percentile (",
    number(shared_p99, accuracy = 0.01),
    "). Genes beyond the displayed range: ",
    outside_note,
    ". Density and box statistics use all significant genes."
  )
)

save_candidate(
  p_zoom,
  "figure3D_raincloud_candidate_2_publication_zoom"
)

p_pseudolog <- make_raincloud(
  data = lfc,
  points = point_data,
  x_trans = "pseudo_log",
  title = "Magnitude of significant transcriptional responses",
  subtitle = "Pseudo-log scale reveals central distributions and extreme responses",
  caption = paste0(
    "Points show a deterministic sample of up to ",
    point_sample_size,
    " genes per fungus; clouds and box summaries use all significant genes."
  )
)

save_candidate(
  p_pseudolog,
  "figure3D_raincloud_candidate_3_pseudolog"
)

message("\nGenerated raincloud candidates:")
message("  1. Full linear range")
message("  2. Publication zoom to shared 99th percentile")
message("  3. Pseudo-log exploratory scale")
message("\nFiles written to:")
message("  ", output_dir)
message("\nDistribution summary:")
print(summary_data, n = Inf)
