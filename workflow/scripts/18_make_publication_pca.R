#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(patchwork)
})

# =============================================================================
# Final publication PCA figure
#
# Two PCA panels:
#   A. Fusarium cf. salinense
#   B. Ophiostoma novo-ulmi
#
# One custom horizontal legend is placed beneath both panels.
# =============================================================================

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

fusarium_coordinates_file <- file.path(
  "results",
  "fusarium",
  "deseq2_results",
  "fusarium_pca_coordinates.tsv"
)

fusarium_variance_file <- file.path(
  "results",
  "fusarium",
  "deseq2_results",
  "fusarium_pca_variance.tsv"
)

ophiostoma_coordinates_file <- file.path(
  "results",
  "ophiostoma",
  "deseq2_qc",
  "pca_coordinates.tsv"
)

ophiostoma_vst_file <- file.path(
  "results",
  "ophiostoma",
  "deseq2_qc",
  "ophiostoma_gene_level_vst_blind.rds"
)

figure_output_dir <- file.path(
  "figures",
  "publication"
)

results_output_dir <- file.path(
  "results",
  "publication",
  "figure3",
  "pca"
)

dir.create(
  figure_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  results_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# -----------------------------------------------------------------------------
# Validate inputs
# -----------------------------------------------------------------------------

required_files <- c(
  fusarium_coordinates_file,
  fusarium_variance_file,
  ophiostoma_coordinates_file,
  ophiostoma_vst_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    paste0(
      "Missing required file(s):\n",
      paste0("  - ", missing_files, collapse = "\n")
    )
  )
}

# -----------------------------------------------------------------------------
# Final trio4 condition mapping
# -----------------------------------------------------------------------------

condition_levels <- c(
  "Interaction zone",
  "Self-interaction control",
  "Outside interaction zone"
)

condition_colours <- c(
  "Interaction zone" = "#94475E",
  "Self-interaction control" = "#364C54",
  "Outside interaction zone" = "#E5A11F"
)

condition_shapes <- c(
  "Interaction zone" = 21,
  "Self-interaction control" = 24,
  "Outside interaction zone" = 22
)

# -----------------------------------------------------------------------------
# Read Fusarium PCA
# -----------------------------------------------------------------------------

fusarium_pca <- read.delim(
  fusarium_coordinates_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

fusarium_variance <- read.delim(
  fusarium_variance_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

fusarium_pc1 <- fusarium_variance$percent_variance[
  fusarium_variance$component == "PC1"
]

fusarium_pc2 <- fusarium_variance$percent_variance[
  fusarium_variance$component == "PC2"
]

if (
  length(fusarium_pc1) != 1L ||
  length(fusarium_pc2) != 1L
) {
  stop("Could not recover Fusarium PC1 and PC2 variance values.")
}

fusarium_label_map <- c(
  interaction = "Interaction zone",
  control = "Self-interaction control"
)

unknown_fusarium_conditions <- setdiff(
  unique(fusarium_pca$condition),
  names(fusarium_label_map)
)

if (length(unknown_fusarium_conditions) > 0L) {
  stop(
    paste0(
      "Unknown Fusarium condition(s): ",
      paste(unknown_fusarium_conditions, collapse = ", ")
    )
  )
}

fusarium_pca$condition_label <- unname(
  fusarium_label_map[fusarium_pca$condition]
)

fusarium_pca$condition_label <- factor(
  fusarium_pca$condition_label,
  levels = condition_levels
)

# -----------------------------------------------------------------------------
# Read Ophiostoma PCA
# -----------------------------------------------------------------------------

ophiostoma_pca <- read.delim(
  ophiostoma_coordinates_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ophiostoma_label_map <- c(
  interaction = "Interaction zone",
  self = "Outside interaction zone",
  onu = "Self-interaction control"
)

unknown_ophiostoma_conditions <- setdiff(
  unique(ophiostoma_pca$condition),
  names(ophiostoma_label_map)
)

if (length(unknown_ophiostoma_conditions) > 0L) {
  stop(
    paste0(
      "Unknown Ophiostoma condition(s): ",
      paste(unknown_ophiostoma_conditions, collapse = ", ")
    )
  )
}

ophiostoma_pca$condition_label <- unname(
  ophiostoma_label_map[ophiostoma_pca$condition]
)

ophiostoma_pca$condition_label <- factor(
  ophiostoma_pca$condition_label,
  levels = condition_levels
)

# -----------------------------------------------------------------------------
# Recover Ophiostoma PCA variance
# -----------------------------------------------------------------------------

ophiostoma_vst <- readRDS(
  ophiostoma_vst_file
)

ophiostoma_plotpca_data <- DESeq2::plotPCA(
  ophiostoma_vst,
  intgroup = "condition",
  returnData = TRUE
)

ophiostoma_percent_var <- round(
  100 * attr(
    ophiostoma_plotpca_data,
    "percentVar"
  ),
  1
)

ophiostoma_pc1 <- ophiostoma_percent_var[1]
ophiostoma_pc2 <- ophiostoma_percent_var[2]

# -----------------------------------------------------------------------------
# Publication theme
# -----------------------------------------------------------------------------

theme_publication_pca <- function() {
  theme_classic(base_size = 10, base_family = "Nimbus Sans") +
    theme(
      axis.title = element_text(
        colour = "#222222",
        size = 10
      ),
      axis.text = element_text(
        colour = "#333333",
        size = 8.8
      ),
      axis.line = element_line(
        colour = "#333333",
        linewidth = 0.45
      ),
      axis.ticks = element_line(
        colour = "#333333",
        linewidth = 0.45
      ),
      plot.title = element_text(
        size = 11,
        face = "italic",
        hjust = 0.5,
        margin = margin(b = 7)
      ),
      plot.tag = element_text(
        face = "bold",
        size = 10,
        colour = "#222222"
      ),
      plot.tag.position = c(0.015, 1.02),
      legend.position = "none",
      plot.margin = margin(
        9,
        10,
        8,
        9
      )
    )
}

# -----------------------------------------------------------------------------
# PCA plot function
# -----------------------------------------------------------------------------

make_pca_plot <- function(
    data,
    pc1,
    pc2,
    species_title,
    panel_tag) {

  ggplot(
    data = data,
    mapping = aes(
      x = PC1,
      y = PC2,
      fill = condition_label,
      shape = condition_label
    )
  ) +
    geom_hline(
      yintercept = 0,
      colour = "#D9D9D9",
      linewidth = 0.45
    ) +
    geom_vline(
      xintercept = 0,
      colour = "#D9D9D9",
      linewidth = 0.45
    ) +
    geom_point(
      size = 5,
      stroke = 0.45,
      colour = "#333333"
    ) +
    scale_fill_manual(
      values = condition_colours,
      limits = condition_levels,
      drop = FALSE
    ) +
    scale_shape_manual(
      values = condition_shapes,
      limits = condition_levels,
      drop = FALSE
    ) +
    labs(
      x = paste0(
        "PC1 (",
        pc1,
        "%)"
      ),
      y = paste0(
        "PC2 (",
        pc2,
        "%)"
      ),
      tag = panel_tag
    ) +
    coord_cartesian(
      clip = "off"
    ) +
    theme_publication_pca()
}

# -----------------------------------------------------------------------------
# Build PCA panels
# -----------------------------------------------------------------------------

fusarium_plot <- make_pca_plot(
  data = fusarium_pca,
  pc1 = fusarium_pc1,
  pc2 = fusarium_pc2,
  species_title = "Fusarium cf. salinense",
  panel_tag = "a"
)

ophiostoma_plot <- make_pca_plot(
  data = ophiostoma_pca,
  pc1 = ophiostoma_pc1,
  pc2 = ophiostoma_pc2,
  species_title = "Ophiostoma novo-ulmi",
  panel_tag = "b"
)

panel_row <- (
  fusarium_plot |
    ophiostoma_plot
) +
  plot_layout(
    widths = c(1, 1.10)
  )

# -----------------------------------------------------------------------------
# Build one custom horizontal legend
# -----------------------------------------------------------------------------

legend_data <- data.frame(
  condition_label = factor(
    condition_levels,
    levels = condition_levels
  ),
  x = c(1, 2, 3),
  y = 1
)

legend_panel <- ggplot(
  legend_data,
  aes(
    x = x,
    y = y,
    fill = condition_label,
    shape = condition_label
  )
) +
  geom_point(
    size = 4,
    stroke = 0.45,
    colour = "#333333"
  ) +
  geom_text(
    aes(label = condition_label),
    nudge_x = 0.11,
    hjust = 0,
    size = 3.05,
    family = "Nimbus Sans",
    colour = "#222222"
  ) +
  scale_fill_manual(
    values = condition_colours,
    limits = condition_levels,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = condition_shapes,
    limits = condition_levels,
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(0.72, 4.02),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0.78, 1.22),
    expand = c(0, 0)
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  theme_void(base_family = "Nimbus Sans") +
  theme(
    legend.position = "none",
    plot.margin = margin(
      -2,
      4,
      0,
      4
    )
  )

combined_plot <- (
  panel_row /
    legend_panel
) +
  plot_layout(
    heights = c(1, 0.075)
  )

# -----------------------------------------------------------------------------
# Save final outputs
# -----------------------------------------------------------------------------

pdf_file <- file.path(
  figure_output_dir,
  "figure_3ab_pca.pdf"
)

png_file <- file.path(
  figure_output_dir,
  "figure_3ab_pca.png"
)

ggsave(
  filename = pdf_file,
  plot = combined_plot,
  width = 170,
  height = 87.1,
  units = "mm",
  device = grDevices::cairo_pdf,
  bg = "white"
)

ggsave(
  filename = png_file,
  plot = combined_plot,
  width = 170,
  height = 87.1,
  units = "mm",
  dpi = 600,
  bg = "white"
)

# -----------------------------------------------------------------------------
# Save mapping for manuscript records
# -----------------------------------------------------------------------------

palette_table <- data.frame(
  condition = condition_levels,
  shape = unname(
    condition_shapes[condition_levels]
  ),
  hex = unname(
    condition_colours[condition_levels]
  ),
  stringsAsFactors = FALSE
)

write.table(
  palette_table,
  file = file.path(
    results_output_dir,
    "figure_3ab_pca_palette.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# -----------------------------------------------------------------------------
# Save PCA source data for manuscript records
# -----------------------------------------------------------------------------

write.table(
  fusarium_pca,
  file = file.path(
    results_output_dir,
    "fusarium_pca_coordinates.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    component = c("PC1", "PC2"),
    percent_variance = c(
      fusarium_pc1,
      fusarium_pc2
    )
  ),
  file = file.path(
    results_output_dir,
    "fusarium_pca_variance.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  ophiostoma_pca,
  file = file.path(
    results_output_dir,
    "ophiostoma_pca_coordinates.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    component = c("PC1", "PC2"),
    percent_variance = c(
      ophiostoma_pc1,
      ophiostoma_pc2
    )
  ),
  file = file.path(
    results_output_dir,
    "ophiostoma_pca_variance.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Wrote ", pdf_file)
message("Wrote ", png_file)
message(
  "Fusarium variance: PC1 = ",
  fusarium_pc1,
  "%; PC2 = ",
  fusarium_pc2,
  "%."
)
message(
  "Ophiostoma variance: PC1 = ",
  ophiostoma_pc1,
  "%; PC2 = ",
  ophiostoma_pc2,
  "%."
)
message("Final PCA figure completed successfully.")
