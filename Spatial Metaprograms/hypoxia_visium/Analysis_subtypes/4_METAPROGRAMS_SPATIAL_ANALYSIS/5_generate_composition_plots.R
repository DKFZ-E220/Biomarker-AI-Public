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
# COMPOSITION PLOTS: Generate 3 comprehensive composition visualizations
# ==============================================================================
# Output:
# 1. CellLine_MP_Composition_Radioresistance.pdf - Cell line stacked bars ordered by radioresistance
# 2. Sample_MP_Composition_Radioresistance.pdf - All samples ordered by cell line radioresistance
# 3. All_Samples_Hypoxia_MP_Composition.pdf - Large format integrated view
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Wong colorblind-safe palette
wong_palette <- list(
  MP1 = "#0072B2",      # Blue
  MP2 = "#E69F00",      # Orange
  MP3 = "#009E73",      # Green
  MP4 = "#CC79A7",      # Mauve
  MP5 = "#D55E00",      # Vermillion
  Hypoxic = "#E74C3C",  # Red
  Normoxic = "#3498DB"  # Light Blue
)

# Radioresistance ranking (TC50 order: least to most sensitive)
cell_line_ranking <- c("FaDu", "SAS", "UT5", "UT8", "Cal33", "SAT")

# Data paths
BASE_DIR <- file.path(dirname(dirname(analysis_subtypes_dir)))
MASTER_SPOT_TABLE <- file.path(analysis_subtypes_dir, "2_METAPROGRAMS_INTEGRATION",
                                "data", "master_spot_table.csv")

# Output directory
OUTPUT_PARENT_DIR <- file.path(BASE_DIR, "Writing")
dir.create(OUTPUT_PARENT_DIR, showWarnings = FALSE, recursive = TRUE)

# Publication theme
theme_publication <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
}

cat("═══════════════════════════════════════════════════════════════\n")
cat("   COMPOSITION PLOTS: METAPROGRAM COMPOSITION BY CELL LINE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading data...\n")

if (!file.exists(MASTER_SPOT_TABLE)) {
  stop("Master spot table not found:", MASTER_SPOT_TABLE)
}

master_data <- read_csv(MASTER_SPOT_TABLE, show_col_types = FALSE)
cat("Master data:", nrow(master_data), "spots\n")
cat("Cell lines:", n_distinct(master_data$cell_line), "\n\n")

# Prepare data
composition_data <- master_data %>%
  mutate(
    hypoxia_status = case_when(
      hypoxia_label == 1  ~ "Hypoxic",
      hypoxia_label == 0  ~ "Normoxic",
      TRUE                ~ NA_character_
    ),
    cell_line = factor(cell_line, levels = cell_line_ranking)
  ) %>%
  filter(!is.na(hypoxia_status), !is.na(dominant_metaprogram)) %>%
  as.data.frame()

mp_colors <- c(
  "MP1" = wong_palette$MP1,
  "MP2" = wong_palette$MP2,
  "MP3" = wong_palette$MP3,
  "MP4" = wong_palette$MP4,
  "MP5" = wong_palette$MP5
)

# ==============================================================================
# PLOT 1: Cell Line Composition (Ordered by Radioresistance)
# ==============================================================================

cat("Plot 1: Cell Line MP Composition (Radioresistance Ordered)...\n")

cellline_composition <- composition_data %>%
  group_by(cell_line, dominant_metaprogram, hypoxia_status) %>%
  summarise(n_spots = n(), .groups = "drop") %>%
  group_by(cell_line, hypoxia_status) %>%
  mutate(percentage = n_spots / sum(n_spots) * 100) %>%
  ungroup()

p1 <- ggplot(cellline_composition, aes(x = factor(cell_line, levels = cell_line_ranking), 
                                        y = percentage, 
                                        fill = dominant_metaprogram)) +
  geom_col(alpha = 0.8, color = "black", size = 0.3) +
  facet_wrap(~hypoxia_status, ncol = 2) +
  scale_fill_manual(
    values = mp_colors,
    name = "Dominant MP"
  ) +
  labs(
    title = "Metaprogram Composition by Cell Line (Ordered by Radioresistance)",
    subtitle = "Left: Normoxic | Right: Hypoxic",
    x = "Cell Line (← more radioresistant | more radiosensitive →)",
    y = "Percentage of Spots (%)"
  ) +
  theme_publication() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 11, face = "bold")
  )

ggsave(
  file.path(OUTPUT_PARENT_DIR, "CellLine_MP_Composition_Radioresistance.pdf"),
  p1,
  width = 12,
  height = 7,
  dpi = 300
)

cat("  Saved: CellLine_MP_Composition_Radioresistance.pdf\n")

# ==============================================================================
# PLOT 2: Sample Composition (Ordered by Cell Line Radioresistance)
# ==============================================================================

cat("Plot 2: Sample MP Composition (Radioresistance Ordered)...\n")

sample_composition <- composition_data %>%
  group_by(sample, cell_line, dominant_metaprogram, hypoxia_status) %>%
  summarise(n_spots = n(), .groups = "drop") %>%
  group_by(sample, cell_line, hypoxia_status) %>%
  mutate(percentage = n_spots / sum(n_spots) * 100) %>%
  ungroup() %>%
  mutate(
    cell_line = factor(cell_line, levels = cell_line_ranking),
    sample = reorder(sample, as.numeric(cell_line))
  )

p2 <- ggplot(sample_composition, aes(x = sample, y = percentage, fill = dominant_metaprogram)) +
  geom_col(alpha = 0.8, color = "black", size = 0.3) +
  facet_wrap(~hypoxia_status, ncol = 2) +
  scale_fill_manual(
    values = mp_colors,
    name = "Dominant MP"
  ) +
  labs(
    title = "Metaprogram Composition by Sample (Cell Lines Ordered by Radioresistance)",
    subtitle = "Left: Normoxic | Right: Hypoxic",
    x = "Sample",
    y = "Percentage of Spots (%)"
  ) +
  theme_publication() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    strip.text = element_text(size = 11, face = "bold")
  )

ggsave(
  file.path(OUTPUT_PARENT_DIR, "Sample_MP_Composition_Radioresistance.pdf"),
  p2,
  width = 14,
  height = 7,
  dpi = 300
)

cat("  Saved: Sample_MP_Composition_Radioresistance.pdf\n")

# ==============================================================================
# PLOT 3: All Samples Integrated Hypoxia-MP Composition
# ==============================================================================

cat("Plot 3: All Samples Integrated Hypoxia MP Composition...\n")

all_samples_composition <- composition_data %>%
  group_by(sample, dominant_metaprogram, hypoxia_status) %>%
  summarise(n_spots = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = hypoxia_status,
    values_from = n_spots,
    values_fill = 0
  )

all_samples_long <- all_samples_composition %>%
  pivot_longer(
    cols = c("Normoxic", "Hypoxic"),
    names_to = "hypoxia_status",
    values_to = "n_spots"
  ) %>%
  mutate(
    sample = factor(sample),
    hypoxia_status = factor(hypoxia_status, levels = c("Normoxic", "Hypoxic"))
  )

p3 <- ggplot(all_samples_long, aes(x = sample, y = n_spots, fill = dominant_metaprogram)) +
  geom_col(alpha = 0.8, color = "black", size = 0.2) +
  facet_wrap(~hypoxia_status, ncol = 2, scales = "free_y") +
  scale_fill_manual(
    values = mp_colors,
    name = "Dominant MP"
  ) +
  labs(
    title = "Metaprogram Composition: All Samples Integrated Analysis",
    subtitle = "Left: Normoxic Regions | Right: Hypoxic Regions",
    x = "Sample",
    y = "Number of Spots"
  ) +
  theme_publication() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
    strip.text = element_text(size = 11, face = "bold")
  )

ggsave(
  file.path(OUTPUT_PARENT_DIR, "All_Samples_Hypoxia_MP_Composition.pdf"),
  p3,
  width = 16,
  height = 9,
  dpi = 300
)

cat("  Saved: All_Samples_Hypoxia_MP_Composition.pdf\n\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("COMPOSITION PLOTS COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Generated 3 composition plots:\n")
cat("  1. CellLine_MP_Composition_Radioresistance.pdf\n")
cat("  2. Sample_MP_Composition_Radioresistance.pdf\n")
cat("  3. All_Samples_Hypoxia_MP_Composition.pdf\n\n")

cat("Output directory:", OUTPUT_PARENT_DIR, "\n\n")

cat("Cell line radioresistance order (← more resistant | more sensitive →):\n")
for (i in seq_along(cell_line_ranking)) {
  cat("  ", i, ". ", cell_line_ranking[i], "\n", sep = "")
}

cat("\nWong colorblind-safe palette used:\n")
cat("  MP1: #0072B2 (Blue)\n")
cat("  MP2: #E69F00 (Orange)\n")
cat("  MP3: #009E73 (Green)\n")
cat("  MP4: #CC79A7 (Mauve)\n\n")

