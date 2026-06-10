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
# QC SUMMARY PLOTS: Generate comprehensive QC visualizations for all cell lines
# ==============================================================================
# Generates 6 PDF plots with Wong colorblind-safe palette:
# 1. Hypoxia Status Distribution (counts)
# 2. MP Score Distributions (boxplots)
# 3. Dominant MP Distribution (stacked bar)
# 4. Sample Breakdown (bar plot)
# 5. Spot Count Summary (bar plot)
# 6. Hypoxia Percentage Distribution (violin + box)
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(gridExtra)

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

# Data paths
BASE_DIR <- file.path(dirname(dirname(analysis_subtypes_dir)))
MASTER_SPOT_TABLE <- file.path(analysis_subtypes_dir, "2_METAPROGRAMS_INTEGRATION", 
                                "data", "master_spot_table.csv")

# Output directory
OUTPUT_PARENT_DIR <- file.path(BASE_DIR, "Writing")
OUTPUT_DIR <- file.path(OUTPUT_PARENT_DIR, "QC_Summary_Plots_AllCellLines")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

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
      panel.grid.major = element_line(color = "gray95", size = 0.2),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
}

cat("═══════════════════════════════════════════════════════════════\n")
cat("   QC SUMMARY PLOTS: ALL CELL LINES COMBINED\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Output directory:", OUTPUT_DIR, "\n\n")

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading data...\n")

if (!file.exists(MASTER_SPOT_TABLE)) {
  stop("Master spot table not found:", MASTER_SPOT_TABLE)
}

master_data <- read_csv(MASTER_SPOT_TABLE, show_col_types = FALSE)
cat("Master data:", nrow(master_data), "spots\n")
cat("Cell lines:", n_distinct(master_data$cell_line), "\n")
cat("Samples:", n_distinct(master_data$sample), "\n\n")

# Prepare data
qc_data <- master_data %>%
  mutate(
    hypoxia_status = case_when(
      hypoxia_label == 1  ~ "Hypoxic",
      hypoxia_label == 0  ~ "Normoxic",
      TRUE                ~ NA_character_
    ),
    cell_line = factor(cell_line),
    sample = factor(sample)
  ) %>%
  filter(!is.na(hypoxia_status), !is.na(dominant_metaprogram)) %>%
  as.data.frame()

# ==============================================================================
# PLOT 1: Hypoxia Status Distribution
# ==============================================================================

cat("Plot 1: Hypoxia Status Distribution...\n")

hypoxia_counts <- qc_data %>%
  group_by(cell_line, hypoxia_status) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = count / sum(count) * 100)

p1 <- ggplot(hypoxia_counts, aes(x = cell_line, y = count, fill = hypoxia_status)) +
  geom_col(position = "dodge", alpha = 0.8, color = "black", size = 0.3) +
  scale_fill_manual(
    values = c("Hypoxic" = wong_palette$Hypoxic, "Normoxic" = wong_palette$Normoxic),
    name = "Status"
  ) +
  labs(
    title = "Hypoxia Status Distribution Across Cell Lines",
    x = "Cell Line",
    y = "Number of Spots"
  ) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(OUTPUT_DIR, "01_Hypoxia_Status_Distribution.pdf"),
  p1,
  width = 10,
  height = 6,
  dpi = 300
)

cat("  Saved: 01_Hypoxia_Status_Distribution.pdf\n")

# ==============================================================================
# PLOT 2: MP Score Distributions (Boxplots)
# ==============================================================================

cat("Plot 2: MP Score Distributions...\n")

mp_long <- qc_data %>%
  pivot_longer(
    cols = c("MP1", "MP2", "MP3", "MP4", "MP5"),
    names_to = "Metaprogram",
    values_to = "Score"
  )

mp_order <- c("MP1", "MP2", "MP3", "MP4", "MP5")
mp_colors <- c(
  "MP1" = wong_palette$MP1,
  "MP2" = wong_palette$MP2,
  "MP3" = wong_palette$MP3,
  "MP4" = wong_palette$MP4,
  "MP5" = wong_palette$MP5
)

p2 <- ggplot(mp_long, aes(x = factor(Metaprogram, levels = mp_order), y = Score, fill = Metaprogram)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1, color = "black", size = 0.3) +
  scale_fill_manual(values = mp_colors, name = "Metaprogram", guide = "none") +
  labs(
    title = "Metaprogram Score Distributions Across All Spots",
    x = "Metaprogram",
    y = "MP Score"
  ) +
  theme_publication()

ggsave(
  file.path(OUTPUT_DIR, "02_MP_Score_Distributions.pdf"),
  p2,
  width = 10,
  height = 6,
  dpi = 300
)

cat("  Saved: 02_MP_Score_Distributions.pdf\n")

# ==============================================================================
# PLOT 3: Dominant MP Distribution (Stacked Bar)
# ==============================================================================

cat("Plot 3: Dominant MP Distribution...\n")

dominant_counts <- qc_data %>%
  group_by(cell_line, dominant_metaprogram) %>%
  summarise(count = n(), .groups = "drop")

p3 <- ggplot(dominant_counts, aes(x = cell_line, y = count, fill = dominant_metaprogram)) +
  geom_col(alpha = 0.8, color = "black", size = 0.3) +
  scale_fill_manual(
    values = mp_colors,
    name = "Dominant MP"
  ) +
  labs(
    title = "Dominant Metaprogram Distribution by Cell Line",
    x = "Cell Line",
    y = "Number of Spots"
  ) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(OUTPUT_DIR, "03_Dominant_MP_Distribution.pdf"),
  p3,
  width = 10,
  height = 6,
  dpi = 300
)

cat("  Saved: 03_Dominant_MP_Distribution.pdf\n")

# ==============================================================================
# PLOT 4: Sample Breakdown (Samples per Cell Line)
# ==============================================================================

cat("Plot 4: Sample Breakdown...\n")

sample_counts <- qc_data %>%
  group_by(cell_line, sample) %>%
  summarise(n_spots = n(), .groups = "drop") %>%
  group_by(cell_line) %>%
  summarise(n_samples = n_distinct(sample), .groups = "drop")

p4 <- ggplot(sample_counts, aes(x = reorder(cell_line, n_samples), y = n_samples)) +
  geom_col(fill = wong_palette$MP1, alpha = 0.8, color = "black", size = 0.3) +
  labs(
    title = "Number of Samples per Cell Line",
    x = "Cell Line",
    y = "Number of Samples"
  ) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(OUTPUT_DIR, "04_Sample_Breakdown.pdf"),
  p4,
  width = 10,
  height = 6,
  dpi = 300
)

cat("  Saved: 04_Sample_Breakdown.pdf\n")

# ==============================================================================
# PLOT 5: Spot Count Summary
# ==============================================================================

cat("Plot 5: Spot Count Summary...\n")

spot_counts <- qc_data %>%
  group_by(cell_line) %>%
  summarise(n_spots = n(), .groups = "drop")

p5 <- ggplot(spot_counts, aes(x = reorder(cell_line, n_spots), y = n_spots)) +
  geom_col(fill = wong_palette$MP2, alpha = 0.8, color = "black", size = 0.3) +
  labs(
    title = "Total Spot Counts per Cell Line",
    x = "Cell Line",
    y = "Number of Spots"
  ) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(OUTPUT_DIR, "05_Spot_Count_Summary.pdf"),
  p5,
  width = 10,
  height = 6,
  dpi = 300
)

cat("  Saved: 05_Spot_Count_Summary.pdf\n")

# ==============================================================================
# PLOT 6: Hypoxia Percentage Distribution (Violin + Box)
# ==============================================================================

cat("Plot 6: Hypoxia Percentage Distribution...\n")

hypoxia_pct <- qc_data %>%
  group_by(cell_line, sample) %>%
  summarise(
    hypoxia_pct = sum(hypoxia_label == 1) / n() * 100,
    .groups = "drop"
  )

p6 <- ggplot(hypoxia_pct, aes(x = reorder(cell_line, hypoxia_pct, FUN = median), y = hypoxia_pct)) +
  geom_violin(fill = wong_palette$Hypoxic, alpha = 0.4, color = wong_palette$Hypoxic) +
  geom_boxplot(width = 0.2, fill = wong_palette$Hypoxic, alpha = 0.8, color = "black", size = 0.3) +
  geom_point(alpha = 0.5, size = 2, color = wong_palette$Hypoxic) +
  labs(
    title = "Hypoxia Percentage Distribution per Sample",
    x = "Cell Line",
    y = "Hypoxic Spots (%)"
  ) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(OUTPUT_DIR, "06_Hypoxia_Percentage_Distribution.pdf"),
  p6,
  width = 10,
  height = 6,
  dpi = 300
)

cat("  Saved: 06_Hypoxia_Percentage_Distribution.pdf\n\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("QC SUMMARY PLOTS COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Generated 6 QC plots:\n")
cat("  1. 01_Hypoxia_Status_Distribution.pdf\n")
cat("  2. 02_MP_Score_Distributions.pdf\n")
cat("  3. 03_Dominant_MP_Distribution.pdf\n")
cat("  4. 04_Sample_Breakdown.pdf\n")
cat("  5. 05_Spot_Count_Summary.pdf\n")
cat("  6. 06_Hypoxia_Percentage_Distribution.pdf\n\n")

cat("Output directory:", OUTPUT_DIR, "\n\n")

cat("Wong colorblind-safe palette used:\n")
cat("  MP1: #0072B2 (Blue)\n")
cat("  MP2: #E69F00 (Orange)\n")
cat("  MP3: #009E73 (Green)\n")
cat("  MP4: #CC79A7 (Mauve)\n")
cat("  Hypoxic: #E74C3C (Red)\n")
cat("  Normoxic: #3498DB (Light Blue)\n\n")

