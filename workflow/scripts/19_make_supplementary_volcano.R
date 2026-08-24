#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# =============================================================================
# Supplementary volcano plots
#
# Panel A:
#   Fusarium interaction zone vs self-interaction control
#   x-axis uses apeglm-shrunken log2 fold changes.
#
# Panels B-D:
#   Ophiostoma contrasts
#   x-axis uses canonical DESeq2 log2 fold changes.
#
# Significance thresholds:
#   adjusted p-value < 0.05
#   strong differential expression: |log2 fold change| > 1
# =============================================================================

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

fusarium_file <- file.path(
  "results",
  "fusarium",
  "deseq2_results",
  "fusarium_interaction_vs_self_lfcshrunk.tsv"
)

ophiostoma_files <- c(
  interaction_vs_self = file.path(
    "results",
    "ophiostoma",
    "deseq2_results",
    "tables",
    "interaction_vs_self_all_genes.tsv"
  ),
  onu_vs_self = file.path(
    "results",
    "ophiostoma",
    "deseq2_results",
    "tables",
    "onu_vs_self_all_genes.tsv"
  ),
  interaction_vs_onu = file.path(
    "results",
    "ophiostoma",
    "deseq2_results",
    "tables",
    "interaction_vs_onu_all_genes.tsv"
  )
)

figure_output_dir <- file.path(
  "figures",
  "publication"
)

results_output_dir <- file.path(
  "results",
  "publication",
  "supplementary",
  "figure_s01_volcano"
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

required_files <- c(
  fusarium_file,
  unname(ophiostoma_files)
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
# Plotting constants
# -----------------------------------------------------------------------------

alpha_threshold <- 0.05
lfc_threshold <- 1

# remains palette from ltc:
# #78477D  #F0EFC6  #FF7C35  #F0D866
fusarium_accent <- "#FF7C35"
ophiostoma_accent <- "#78477D"

# Neutral background keeps the volcano plots quiet and supplementary.
colour_background <- "#D9D9D9"

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

validate_result_table <- function(data, table_name) {

  required_columns <- c(
    "gene_id",
    "log2FoldChange",
    "padj"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(data)
  )

  if (length(missing_columns) > 0L) {
    stop(
      paste0(
        table_name,
        " is missing required column(s): ",
        paste(missing_columns, collapse = ", ")
      )
    )
  }

  invisible(TRUE)
}

prepare_volcano_data <- function(
    data,
    organism,
    contrast_id,
    lfc_type) {

  validate_result_table(
    data,
    contrast_id
  )

  data$organism <- organism
  data$contrast_id <- contrast_id
  data$lfc_type <- lfc_type

  data$status <- "Background"

  strong <- (
    !is.na(data$padj) &
      data$padj < alpha_threshold &
      !is.na(data$log2FoldChange) &
      abs(data$log2FoldChange) > lfc_threshold
  )

  data$status[strong] <- "Strong differential expression"

  data$status <- factor(
    data$status,
    levels = c(
      "Background",
      "Strong differential expression"
    )
  )

  # Avoid Inf values where padj is exactly zero.
  positive_padj <- data$padj[
    !is.na(data$padj) &
      data$padj > 0
  ]

  if (length(positive_padj) == 0L) {
    stop(
      "No positive adjusted p-values were found for ",
      contrast_id,
      "."
    )
  }

  minimum_positive_padj <- min(
    positive_padj
  )

  data$padj_for_plot <- data$padj

  data$padj_for_plot[
    !is.na(data$padj_for_plot) &
      data$padj_for_plot <= 0
  ] <- minimum_positive_padj

  data$minus_log10_padj <- -log10(
    data$padj_for_plot
  )

  data
}

make_summary_row <- function(data) {

  data.frame(
    organism = unique(data$organism),
    contrast = unique(data$contrast_id),
    lfc_type = unique(data$lfc_type),
    total_rows = nrow(data),
    genes_with_padj = sum(!is.na(data$padj)),
    significant_padj_lt_0_05 = sum(
      !is.na(data$padj) &
        data$padj < alpha_threshold
    ),
    strong_padj_lt_0_05_abs_lfc_gt_1 = sum(
      !is.na(data$padj) &
        data$padj < alpha_threshold &
        !is.na(data$log2FoldChange) &
        abs(data$log2FoldChange) > lfc_threshold
    ),
    strong_up = sum(
      !is.na(data$padj) &
        data$padj < alpha_threshold &
        data$log2FoldChange > lfc_threshold
    ),
    strong_down = sum(
      !is.na(data$padj) &
        data$padj < alpha_threshold &
        data$log2FoldChange < -lfc_threshold
    ),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Read and prepare Fusarium data
# -----------------------------------------------------------------------------

fusarium <- read.delim(
  fusarium_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

fusarium <- prepare_volcano_data(
  data = fusarium,
  organism = "Fusarium",
  contrast_id = "interaction_vs_control",
  lfc_type = "apeglm-shrunken"
)

# -----------------------------------------------------------------------------
# Read and prepare Ophiostoma data
# -----------------------------------------------------------------------------

ophiostoma_interaction_self <- read.delim(
  ophiostoma_files[["interaction_vs_self"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ophiostoma_interaction_self <- prepare_volcano_data(
  data = ophiostoma_interaction_self,
  organism = "Ophiostoma",
  contrast_id = "interaction_vs_self",
  lfc_type = "unshrunk DESeq2"
)

ophiostoma_onu_self <- read.delim(
  ophiostoma_files[["onu_vs_self"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ophiostoma_onu_self <- prepare_volcano_data(
  data = ophiostoma_onu_self,
  organism = "Ophiostoma",
  contrast_id = "onu_vs_self",
  lfc_type = "unshrunk DESeq2"
)

ophiostoma_interaction_onu <- read.delim(
  ophiostoma_files[["interaction_vs_onu"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ophiostoma_interaction_onu <- prepare_volcano_data(
  data = ophiostoma_interaction_onu,
  organism = "Ophiostoma",
  contrast_id = "interaction_vs_onu",
  lfc_type = "unshrunk DESeq2"
)

all_datasets <- list(
  fusarium = fusarium,
  ophiostoma_interaction_self = ophiostoma_interaction_self,
  ophiostoma_onu_self = ophiostoma_onu_self,
  ophiostoma_interaction_onu = ophiostoma_interaction_onu
)

# -----------------------------------------------------------------------------
# Shared axis limits
#
# Symmetric x-limits are calculated from all four panels.
# The y-axis is capped at the 99.9th percentile so only the most extreme
# adjusted p-values are clipped, reducing the visible horizontal ceiling.
# Values above the cap are drawn at the upper boundary.
# -----------------------------------------------------------------------------

all_absolute_lfc <- unlist(
  lapply(
    all_datasets,
    function(data) {
      abs(
        data$log2FoldChange[
          is.finite(data$log2FoldChange)
        ]
      )
    }
  )
)

x_limit <- ceiling(
  as.numeric(
    quantile(
      all_absolute_lfc,
      probs = 0.999,
      na.rm = TRUE,
      names = FALSE
    )
  )
)

x_limit <- max(
  x_limit,
  2
)

all_y_values <- unlist(
  lapply(
    all_datasets,
    function(data) {
      data$minus_log10_padj[
        is.finite(data$minus_log10_padj)
      ]
    }
  )
)

y_limit <- ceiling(
  as.numeric(
    quantile(
      all_y_values,
      probs = 0.995,
      na.rm = TRUE,
      names = FALSE
    )
  )
)

y_limit <- max(
  y_limit,
  -log10(alpha_threshold) + 1
)

for (dataset_name in names(all_datasets)) {

  all_datasets[[dataset_name]]$x_for_plot <- pmax(
    pmin(
      all_datasets[[dataset_name]]$log2FoldChange,
      x_limit
    ),
    -x_limit
  )

  all_datasets[[dataset_name]]$y_for_plot <- pmin(
    all_datasets[[dataset_name]]$minus_log10_padj,
    y_limit
  )
}

fusarium <- all_datasets$fusarium
ophiostoma_interaction_self <- all_datasets$ophiostoma_interaction_self
ophiostoma_onu_self <- all_datasets$ophiostoma_onu_self
ophiostoma_interaction_onu <- all_datasets$ophiostoma_interaction_onu

# -----------------------------------------------------------------------------
# Theme
# -----------------------------------------------------------------------------

theme_volcano <- function() {

  theme_classic(
    base_size = 9.5
  ) +
    theme(
      axis.title = element_text(
        size = 9.5,
        colour = "#222222"
      ),
      axis.text = element_text(
        size = 8,
        colour = "#333333"
      ),
      axis.line = element_line(
        colour = "#333333",
        linewidth = 0.4
      ),
      axis.ticks = element_line(
        colour = "#333333",
        linewidth = 0.35
      ),
      plot.title = element_blank(),
      plot.subtitle = element_text(
        size = 9,
        hjust = 0,
        colour = "#333333",
        margin = margin(
          b = 7
        )
      ),
      plot.tag = element_text(
        size = 13,
        face = "bold",
        colour = "#222222"
      ),
      plot.tag.position = c(
        0.01,
        1.03
      ),
      legend.position = "none",
      plot.margin = margin(
        9,
        9,
        7,
        8
      )
    )
}

# -----------------------------------------------------------------------------
# Volcano plotting function
# -----------------------------------------------------------------------------

make_volcano_plot <- function(
    data,
    organism_title,
    contrast_subtitle,
    panel_tag,
    accent_colour) {

  plot_data <- data[
    is.finite(data$x_for_plot) &
      is.finite(data$y_for_plot),
    ,
    drop = FALSE
  ]

  plot_colours <- c(
    "Background" = colour_background,
    "Strong differential expression" = accent_colour
  )

  ggplot(
    plot_data,
    aes(
      x = x_for_plot,
      y = y_for_plot,
      colour = status
    )
  ) +
    geom_hline(
      yintercept = -log10(alpha_threshold),
      linetype = "dashed",
      colour = "#777777",
      linewidth = 0.35
    ) +
    geom_vline(
      xintercept = c(
        -lfc_threshold,
        lfc_threshold
      ),
      linetype = "dashed",
      colour = "#777777",
      linewidth = 0.35
    ) +
    geom_point(
      size = 0.7,
      alpha = 0.68,
      stroke = 0
    ) +
    scale_colour_manual(
      values = plot_colours,
      drop = FALSE
    ) +
    scale_x_continuous(
      limits = c(
        -x_limit,
        x_limit
      ),
      expand = expansion(
        mult = c(
          0.02,
          0.02
        )
      )
    ) +
    scale_y_continuous(
      limits = c(
        0,
        y_limit
      ),
      expand = expansion(
        mult = c(
          0,
          0.025
        )
      )
    ) +
    labs(
      x = "log2 fold change",
      y = expression(-log[10]("adjusted p-value")),
      subtitle = contrast_subtitle,
      tag = panel_tag
    ) +
    theme_volcano()
}

# -----------------------------------------------------------------------------
# Build panels
# -----------------------------------------------------------------------------

plot_a <- make_volcano_plot(
  data = fusarium,
  organism_title = NULL,
  contrast_subtitle = "Interaction vs self-interaction control",
  panel_tag = "A",
  accent_colour = fusarium_accent
)

plot_b <- make_volcano_plot(
  data = ophiostoma_interaction_self,
  organism_title = NULL,
  contrast_subtitle = "Interaction vs outside zone",
  panel_tag = "B",
  accent_colour = ophiostoma_accent
)

plot_c <- make_volcano_plot(
  data = ophiostoma_onu_self,
  organism_title = NULL,
  contrast_subtitle = "Self-interaction control vs outside zone",
  panel_tag = "C",
  accent_colour = ophiostoma_accent
)

plot_d <- make_volcano_plot(
  data = ophiostoma_interaction_onu,
  organism_title = NULL,
  contrast_subtitle = "Interaction vs self-interaction control",
  panel_tag = "D",
  accent_colour = ophiostoma_accent
)

panel_grid <- wrap_plots(
  plot_a,
  plot_b,
  plot_c,
  plot_d,
  ncol = 2,
  nrow = 2
)

# -----------------------------------------------------------------------------
# Custom shared legend
# -----------------------------------------------------------------------------

legend_data <- data.frame(
  label = factor(
    c(
      "Background genes",
      "Strong DE — Fusarium",
      "Strong DE — Ophiostoma"
    ),
    levels = c(
      "Background genes",
      "Strong DE — Fusarium",
      "Strong DE — Ophiostoma"
    )
  ),
  x = seq_len(3),
  y = 1
)

legend_colours <- c(
  "Background genes" = colour_background,
  "Strong DE — Fusarium" = fusarium_accent,
  "Strong DE — Ophiostoma" = ophiostoma_accent
)

legend_panel <- ggplot(
  legend_data,
  aes(
    x = x,
    y = y,
    colour = label
  )
) +
  geom_point(
    size = 2.5
  ) +
  geom_text(
    aes(label = label),
    nudge_x = 0.08,
    hjust = 0,
    size = 2.9,
    colour = "#222222"
  ) +
  scale_colour_manual(
    values = legend_colours,
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(
      0.75,
      3.95
    ),
    expand = c(
      0,
      0
    )
  ) +
  scale_y_continuous(
    limits = c(
      0.8,
      1.2
    ),
    expand = c(
      0,
      0
    )
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(
      -3,
      2,
      0,
      2
    )
  )

combined_plot <- wrap_plots(
  wrap_elements(
    full = panel_grid
  ),
  wrap_elements(
    full = legend_panel
  ),
  ncol = 1,
  heights = c(
    1,
    0.065
  )
)

# -----------------------------------------------------------------------------
# Save outputs
# -----------------------------------------------------------------------------

pdf_file <- file.path(
  figure_output_dir,
  "supplementary_figure_s01_volcano.pdf"
)

png_file <- file.path(
  figure_output_dir,
  "supplementary_figure_s01_volcano.png"
)

ggsave(
  filename = pdf_file,
  plot = combined_plot,
  width = 8.3,
  height = 7.0,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

ggsave(
  filename = png_file,
  plot = combined_plot,
  width = 8.3,
  height = 7.0,
  units = "in",
  dpi = 600,
  bg = "white"
)

# -----------------------------------------------------------------------------
# Save diagnostics and summary tables
# -----------------------------------------------------------------------------

summary_table <- do.call(
  rbind,
  lapply(
    list(
      fusarium,
      ophiostoma_interaction_self,
      ophiostoma_onu_self,
      ophiostoma_interaction_onu
    ),
    make_summary_row
  )
)

rownames(summary_table) <- NULL

write.table(
  summary_table,
  file = file.path(
    results_output_dir,
    "supplementary_figure_s01_volcano_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

plot_parameters <- data.frame(
  parameter = c(
    "adjusted_p_value_threshold",
    "absolute_log2fc_threshold",
    "shared_absolute_x_limit",
    "shared_minus_log10_padj_limit",
    "fusarium_log2fc_type",
    "ophiostoma_log2fc_type"
  ),
  value = c(
    alpha_threshold,
    lfc_threshold,
    x_limit,
    y_limit,
    "apeglm-shrunken",
    "unshrunk DESeq2"
  ),
  stringsAsFactors = FALSE
)

write.table(
  plot_parameters,
  file = file.path(
    results_output_dir,
    "supplementary_figure_s01_plot_parameters.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Wrote ", pdf_file)
message("Wrote ", png_file)
message(
  "Shared plotting limits: |log2FC| ≤ ",
  x_limit,
  "; -log10 adjusted p-value ≤ ",
  y_limit,
  "."
)
message("Volcano plot summary:")
print(summary_table)
message("Supplementary volcano figure completed successfully.")
