#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

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

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_pdf <- file.path(
  output_dir,
  "Figure5_CAZyme_class_matrix_candidate.pdf"
)

output_png <- file.path(
  output_dir,
  "Figure5_CAZyme_class_matrix_candidate.png"
)

output_table <- file.path(
  output_dir,
  "Figure5_CAZyme_class_matrix_plotting_table.tsv"
)

class_order <- c(
  "GH",
  "GT",
  "AA",
  "CBM",
  "CE",
  "PL"
)

column_order <- c(
  "Fusarium_up",
  "Fusarium_down",
  "Ophiostoma_up",
  "Ophiostoma_down"
)

column_labels <- c(
  "Fusarium_up" = "\u2191",
  "Fusarium_down" = "\u2193",
  "Ophiostoma_up" = "\u2191",
  "Ophiostoma_down" = "\u2193"
)

if (!file.exists(input_file)) {
  stop("Missing input file: ", input_file)
}

class_data <- read_tsv(
  input_file,
  show_col_types = FALSE,
  progress = FALSE
)

required_columns <- c(
  "class",
  "species",
  "up",
  "down",
  "secreted_up",
  "secreted_down"
)

missing_columns <- setdiff(
  required_columns,
  names(class_data)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

complete_data <- tidyr::expand_grid(
  species = c(
    "Fusarium",
    "Ophiostoma"
  ),
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
        up,
        down,
        secreted_up,
        secreted_down
      ),
      ~ replace_na(.x, 0)
    )
  )

plot_data <- bind_rows(
  complete_data %>%
    transmute(
      species,
      class,
      direction = "up",
      significant = up,
      secreted = secreted_up
    ),
  complete_data %>%
    transmute(
      species,
      class,
      direction = "down",
      significant = down,
      secreted = secreted_down
    )
) %>%
  mutate(
    column = paste(
      species,
      direction,
      sep = "_"
    ),
    column = factor(
      column,
      levels = column_order
    ),
    class = factor(
      class,
      levels = rev(class_order)
    ),
    secreted_fraction = if_else(
      significant > 0,
      secreted / significant,
      NA_real_
    ),
    main_label = if_else(
      significant > 0,
      as.character(significant),
      ""
    ),
    sub_label = if_else(
      significant > 0,
      paste0(secreted),
      ""
    )
  )

write_tsv(
  plot_data %>%
    mutate(
      column = as.character(column),
      class = as.character(class)
    ),
  output_table
)

plot_object <- ggplot(
  plot_data,
  aes(
    x = column,
    y = class
  )
) +
  geom_tile(
    aes(
      fill = secreted_fraction
    ),
    width = 0.88,
    height = 0.88,
    colour = "white",
    linewidth = 1.1
  ) +
  geom_text(
    aes(
      label = main_label
    ),
    size = 5.2,
    fontface = "bold",
    colour = "grey15",
    vjust = 0.72
  ) +
  geom_text(
    aes(
      label = sub_label
    ),
    size = 2.8,
    colour = "grey35",
    vjust = -1.15
  ) +
  annotate(
    "text",
    x = 1.5,
    y = length(class_order) + 0.72,
    label = "A",
    fontface = "bold",
    size = 5.5
  ) +
  annotate(
    "text",
    x = 3.5,
    y = length(class_order) + 0.72,
    label = "B",
    fontface = "bold",
    size = 5.5
  ) +
  geom_vline(
    xintercept = 2.5,
    linewidth = 0.5,
    colour = "grey75"
  ) +
  scale_fill_gradient(
    low = "#F5F2F6",
    high = "#5E3A78",
    limits = c(0, 1),
    na.value = "white",
    name = "Secreted fraction",
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  scale_x_discrete(
    labels = column_labels,
    drop = FALSE
  ) +
  scale_y_discrete(
    drop = FALSE
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(
    base_size = 10,
    base_family = "sans"
  ) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      size = 14,
      colour = "grey20",
      margin = margin(
        t = 5
      )
    ),
    axis.text.y = element_text(
      size = 10,
      colour = "black",
      margin = margin(
        r = 8
      )
    ),
    legend.position = "right",
    legend.title = element_text(
      size = 8
    ),
    legend.text = element_text(
      size = 7
    ),
    panel.grid = element_blank(),
    plot.margin = margin(
      t = 20,
      r = 10,
      b = 8,
      l = 8
    )
  )

ggsave(
  filename = output_pdf,
  plot = plot_object,
  width = 6.4,
  height = 4.7,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = output_png,
  plot = plot_object,
  width = 6.4,
  height = 4.7,
  units = "in",
  dpi = 600,
  bg = "white"
)

cat("\n")
cat("============================================================\n")
cat("CAZYME CLASS TILE MATRIX COMPLETE\n")
cat("============================================================\n\n")

cat("PDF:   ", output_pdf, "\n", sep = "")
cat("PNG:   ", output_png, "\n", sep = "")
cat("Table: ", output_table, "\n", sep = "")
