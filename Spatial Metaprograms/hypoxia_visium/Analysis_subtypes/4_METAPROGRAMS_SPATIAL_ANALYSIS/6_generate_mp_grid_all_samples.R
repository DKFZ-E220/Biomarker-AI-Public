#!/usr/bin/env Rscript

# Resolve script directory — robust for source() and Rscript
script_dir <- tryCatch({
  frames <- sys.frames()
  ofile <- NULL
  for (f in rev(frames)) { if (!is.null(f$ofile)) { ofile <- f$ofile; break } }
  if (!is.null(ofile) && nchar(ofile) > 0) dirname(normalizePath(ofile)) else getwd()
}, error = function(e) getwd())
setwd(script_dir)
# Normalise to Analysis_subtypes regardless of whether script_dir resolved to the R folder or its parent
analysis_subtypes_dir <- if (grepl("4_METAPROGRAMS_SPATIAL_ANALYSIS", getwd(), fixed = TRUE)) dirname(getwd()) else getwd()

# ==============================================================================
# Figure: MP Composition Grid — All Samples (4 columns × 5 rows)
# ==============================================================================
# Shows dominant metaprogram % composition for every sample,
# arranged in a grid ordered by cell line radioresistance.
# Source: master_spot_table.csv
# Output: Writing/All_Samples_MP_Grid.pdf
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------------------

wong_mp <- c(
  MP1 = "#0072B2",   # Blue       — EMT/Mesenchymal
  MP2 = "#E69F00",   # Orange     — Hypoxia Response
  MP3 = "#009E73",   # Green      — Squamous Differentiation
  MP4 = "#CC79A7",   # Mauve      — ECM/Stromal
  MP5 = "#D55E00"    # Vermillion — Keratinocyte/Squamous-epithelial
)

wong_cellline <- c(
  FaDu  = "#0072B2",
  SAS   = "#009E73",
  UT5   = "#D895D0",
  UT8   = "#F0E442",
  Cal33 = "#E69F00",
  SAT   = "#CC79A7"
)

# Radioresistance order (most → least resistant)
cellline_order <- c("FaDu", "SAS", "UT5", "UT8", "Cal33", "SAT")

BASE_DIR <- file.path(dirname(dirname(analysis_subtypes_dir)))
MASTER_SPOT   <- file.path(dirname(getwd()), "2_METAPROGRAMS_INTEGRATION",
                            "data", "master_spot_table.csv")
OUTPUT_DIR    <- file.path(BASE_DIR, "Writing")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------------------------

cat("Loading master_spot_table...\n")
spots <- read_csv(MASTER_SPOT, show_col_types = FALSE)
cat("  Spots:", nrow(spots), "| Samples:", n_distinct(spots$sample), "\n\n")

# ------------------------------------------------------------------------------
# COMPUTE COMPOSITION
# ------------------------------------------------------------------------------

composition <- spots %>%
  count(sample, cell_line, dominant_metaprogram) %>%
  group_by(sample) %>%
  mutate(
    total   = sum(n),
    pct     = n / total * 100,
    cell_line = factor(cell_line, levels = cellline_order)
  ) %>%
  ungroup() %>%
  # Ordered label: cell line + stripped sample ID for readability
  mutate(
    short_id = sub("_SPT.*", "", sub("_[A-Z][a-z0-9]+_", "_", sample)),
    label    = paste0(cell_line, "\n", short_id),
    # Order samples within cell line by total spot count descending
    cell_line_rank = as.integer(cell_line)
  ) %>%
  arrange(cell_line_rank, desc(total))

# Fix factor order for samples: radioresistance order, then within cell line by spot count
sample_levels <- composition %>%
  distinct(sample, label, cell_line_rank, total) %>%
  arrange(cell_line_rank, desc(total)) %>%
  pull(label)

composition <- composition %>%
  mutate(label = factor(label, levels = unique(sample_levels)))

# ------------------------------------------------------------------------------
# PLOT — facet_wrap 4 columns
# ------------------------------------------------------------------------------

theme_grid <- function() {
  theme_minimal(base_size = 9) +
    theme(
      plot.title        = element_text(size = 11, face = "bold", hjust = 0.5),
      strip.text        = element_text(size = 7.5, face = "bold"),
      axis.text.x       = element_blank(),
      axis.ticks.x      = element_blank(),
      axis.title        = element_text(size = 9, face = "bold"),
      legend.title      = element_text(size = 9, face = "bold"),
      legend.text       = element_text(size = 8),
      panel.grid        = element_blank(),
      panel.border      = element_rect(colour = "grey80", fill = NA, linewidth = 0.4)
    )
}

p_grid <- ggplot(composition,
                 aes(x = "", y = pct, fill = dominant_metaprogram)) +
  geom_col(width = 0.8, color = "black", linewidth = 0.25) +
  scale_fill_manual(
    values = wong_mp,
    name   = "Metaprogram",
    labels = c(
      MP1 = "MP1 — EMT/Mesenchymal",
      MP2 = "MP2 — Hypoxia Response",
      MP3 = "MP3 — Squamous Differentiation",
      MP4 = "MP4 — ECM/Stromal",
      MP5 = "MP5 — Keratinocyte/Squamous-epithelial"
    )
  ) +
  facet_wrap(~ label, ncol = 4) +
  labs(
    title = "Metaprogram Composition — All Samples",
    subtitle = "Ordered by radioresistance (FaDu → SAT). Bar = % dominant MP spots.",
    x = NULL,
    y = "% Spots"
  ) +
  theme_grid() +
  guides(fill = guide_legend(keywidth = 0.8, keyheight = 0.8))

ggsave(
  file.path(OUTPUT_DIR, "All_Samples_MP_Grid.pdf"),
  p_grid, width = 10, height = 13, dpi = 300
)
cat("All_Samples_MP_Grid.pdf saved\n")

# ------------------------------------------------------------------------------
# ALSO: one-panel summary bar — all samples side-by-side, ordered by radioresistance
# (horizontal, so sample names readable)
# ------------------------------------------------------------------------------

# Reorder for single-bar plot: use numeric label on x
bar_order <- composition %>%
  distinct(sample, cell_line, cell_line_rank, total) %>%
  arrange(cell_line_rank, desc(total)) %>%
  pull(sample)

composition_bar <- spots %>%
  count(sample, cell_line, dominant_metaprogram) %>%
  group_by(sample) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(
    sample    = factor(sample, levels = bar_order),
    cell_line = factor(cell_line, levels = cellline_order)
  )

# Compute divider positions: after last sample of each cell line (except SAT)
divider_positions <- composition_bar %>%
  distinct(sample, cell_line) %>%
  mutate(sample_rank = as.integer(factor(sample, levels = bar_order))) %>%
  group_by(cell_line) %>%
  summarise(last_rank = max(sample_rank), .groups = "drop") %>%
  filter(cell_line != last(cellline_order[cellline_order %in% unique(composition_bar$cell_line)])) %>%
  pull(last_rank)

p_bar <- ggplot(composition_bar,
                aes(x = sample, y = pct, fill = dominant_metaprogram)) +
  geom_col(color = "black", linewidth = 0.25) +
  scale_fill_manual(
    values = wong_mp,
    name   = "Metaprogram",
    labels = c(
      MP1 = "MP1 — EMT/Mesenchymal",
      MP2 = "MP2 — Hypoxia Response",
      MP3 = "MP3 — Squamous Differentiation",
      MP4 = "MP4 — ECM/Stromal",
      MP5 = "MP5 — Keratinocyte/Squamous-epithelial"
    )
  ) +
  # Vertical dividers between cell lines
  geom_vline(
    xintercept = divider_positions + 0.5,
    linetype = "dashed", color = "grey50", linewidth = 0.4
  ) +
  labs(
    title    = "Metaprogram Composition — All Samples (Radioresistance Ordered)",
    subtitle = "Dashed lines separate cell lines. FaDu (most resistant) → SAT (most sensitive).",
    x        = "Sample",
    y        = "% Spots"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9,  hjust = 0.5, color = "grey40"),
    axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
    axis.title    = element_text(size = 10, face = "bold"),
    legend.title  = element_text(size = 10, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border  = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  ) +
  guides(fill = guide_legend(keywidth = 0.8, keyheight = 0.8))

ggsave(
  file.path(OUTPUT_DIR, "All_Samples_MP_Bar_Ordered.pdf"),
  p_bar, width = 16, height = 6, dpi = 300
)
cat("All_Samples_MP_Bar_Ordered.pdf saved\n")

cat("\nDone. Output:\n")
cat("  Writing/All_Samples_MP_Grid.pdf\n")
cat("  Writing/All_Samples_MP_Bar_Ordered.pdf\n")
