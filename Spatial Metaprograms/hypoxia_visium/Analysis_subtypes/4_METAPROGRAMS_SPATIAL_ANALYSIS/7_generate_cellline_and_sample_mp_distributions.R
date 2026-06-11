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
# Task 3: Generate Cell-Line Comparative Boxplots + Individual Sample MP Plots
# ==============================================================================
# Part A: Comparative boxplots of metaprograms (Hypoxia vs Normoxia) per cell line
# Part B: Individual sample metaprogram distribution plots (all 20 samples)
# Output: One PDF per cell line (Task 3A) + one PDF per sample (Task 3B)
# ==============================================================================

library(tidyverse)
library(Seurat)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(viridis)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Data paths
BASE_DIR <- file.path(dirname(dirname(analysis_subtypes_dir)))
MASTER_SPOT_TABLE <- file.path(analysis_subtypes_dir, "2_METAPROGRAMS_INTEGRATION",
                                "data", "master_spot_table.csv")

# Output directories
OUTPUT_PARENT_DIR <- file.path(BASE_DIR, "Writing")
OUTPUT_3A <- file.path(OUTPUT_PARENT_DIR, "Task3A_CellLine_Comparative_Boxplots")
OUTPUT_3B <- file.path(OUTPUT_PARENT_DIR, "Task3B_Individual_Sample_MP_Distribution")

dir.create(OUTPUT_3A, showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_3B, showWarnings = FALSE, recursive = TRUE)

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
cat("   TASK 3: CELL-LINE COMPARISONS & SAMPLE MP DISTRIBUTIONS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Output directories:\n")
cat("  3A (Comparative):", OUTPUT_3A, "\n")
cat("  3B (Individual) :", OUTPUT_3B, "\n\n")

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading data...\n")

# Load master spot table
if (!file.exists(MASTER_SPOT_TABLE)) {
  stop("Master spot table not found:", MASTER_SPOT_TABLE)
}
master_data <- read_csv(MASTER_SPOT_TABLE, show_col_types = FALSE)
cat("Master data:", nrow(master_data), "spots\n")

# Get unique samples and cell lines
all_samples <- sort(unique(master_data$sample))
all_cell_lines <- sort(unique(master_data$cell_line))
cat("Samples:", length(all_samples), "\n")
cat("Cell lines:", length(all_cell_lines), "\n\n")

# ==============================================================================
# LOAD METAPROGRAMS & COMPUTE SCORES
# ==============================================================================

cat("Setting up metaprograms...\n")

# Metaprograms are already computed in master_spot_table
mp_names <- c("MP1", "MP2", "MP3", "MP4", "MP5")
cat("Metaprograms ready:", paste(mp_names, collapse = ", "), "\n\n")

# Analysis data: metaprograms already computed in master_spot_table
analysis_data <- master_data %>%
  mutate(hypoxia_status = case_when(
    hypoxia_label == 1  ~ "Hypoxic",
    hypoxia_label == 0  ~ "Normoxic",
    TRUE                ~ NA_character_
  )) %>%
  filter(!is.na(hypoxia_status))

cat("Analysis data prepared:", nrow(analysis_data), "spots with MP scores\n\n")

# ==============================================================================
# TASK 3A: CELL-LINE COMPARATIVE BOXPLOTS
# ==============================================================================

cat("TASK 3A: GENERATING CELL-LINE COMPARATIVE BOXPLOTS\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

generated_files_3a <- character()

for (cell_line in all_cell_lines) {
  
  tryCatch({
    cat("  Processing", cell_line, "...")
    
    # Filter data for this cell line
    cell_line_data <- analysis_data %>%
      filter(cell_line == !!cell_line) %>%
      pivot_longer(
        cols = c("MP1", "MP2", "MP3", "MP4", "MP5"),
        names_to = "Metaprogram",
        values_to = "Score"
      )
    
    if (nrow(cell_line_data) == 0) {
      cat(" No data\n")
      next
    }
    
    # Create boxplot: MP × Hypoxia Status
    mp_order <- c("MP1", "MP2", "MP3", "MP4", "MP5")
    mp_labels <- c("MP1\n(EMT/Mesenchymal)", "MP2\n(Hypoxia Response)", "MP3\n(Squamous Diff.)", "MP4\n(ECM/Stromal)", "MP5\n(Keratinocyte/Basal)")
    mp_colors <- c("MP1" = "#0072B2", "MP2" = "#E69F00", "MP3" = "#009E73", "MP4" = "#CC79A7", "MP5" = "#D55E00")

    p_boxplot <- ggplot(cell_line_data, aes(x = Metaprogram, y = Score, fill = hypoxia_status)) +
      geom_boxplot(alpha = 0.8, outlier.size = 1.5, position = position_dodge(width = 0.8)) +
      scale_fill_manual(
        values = c("Hypoxic" = "#E74C3C", "Normoxic" = "#3498DB"),
        name = "Region Type"
      ) +
      scale_x_discrete(
        limits = mp_order,
        labels = mp_labels
      ) +
      labs(
        title = paste("Metaprogram Distribution:", cell_line),
        subtitle = "Hypoxic vs Normoxic Regions",
        x = "Metaprogram",
        y = "MP Score"
      ) +
      theme_publication()
    
    # Summary statistics
    stats_summary <- cell_line_data %>%
      group_by(Metaprogram, hypoxia_status) %>%
      summarise(
        n = n(),
        mean = round(mean(Score, na.rm = TRUE), 3),
        median = round(median(Score, na.rm = TRUE), 3),
        sd = round(sd(Score, na.rm = TRUE), 3),
        .groups = "drop"
      ) %>%
      arrange(Metaprogram, hypoxia_status)
    
    # Create summary table as text
    p_table <- tableGrob(stats_summary, rows = NULL, theme = ttheme_minimal())
    
    # Combine boxplot and table
    combined <- p_boxplot / p_table +
      plot_layout(heights = c(3, 1.5)) +
      plot_annotation(
        theme = theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
      )
    
    # Save PDF
    pdf_filename <- paste0("CellLine_", gsub("\\.", "_", cell_line), "_Comparative_Boxplot.pdf")
    pdf_path <- file.path(OUTPUT_3A, pdf_filename)
    
    ggsave(pdf_path, combined, width = 10, height = 8, dpi = 300)
    cat(" \n")
    
    generated_files_3a <- c(generated_files_3a, pdf_filename)
    
  }, error = function(e) {
    cat(" Error:", e$message, "\n")
  })
}

cat("\nTask 3A complete:", length(generated_files_3a), "cell line plots generated\n\n")

# ==============================================================================
# TASK 3B: INDIVIDUAL SAMPLE MP DISTRIBUTION PLOTS
# ==============================================================================

cat("TASK 3B: GENERATING INDIVIDUAL SAMPLE MP DISTRIBUTION PLOTS\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

generated_files_3b <- character()

for (sample in all_samples) {
  
  tryCatch({
    cat("  Processing", sample, "...")
    
    # Filter data for this sample
    sample_data <- analysis_data %>%
      filter(sample == !!sample) %>%
      pivot_longer(
        cols = c("MP1", "MP2", "MP3", "MP4", "MP5"),
        names_to = "Metaprogram",
        values_to = "Score"
      )
    
    if (nrow(sample_data) == 0) {
      cat(" No data\n")
      next
    }
    
    # Compute mean MP composition
    composition <- sample_data %>%
      group_by(Metaprogram, hypoxia_status) %>%
      summarise(mean_score = mean(Score, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(
        names_from = hypoxia_status,
        values_from = mean_score,
        values_fill = 0
      )
    
    # Create visualizations
    plots <- list()
    
    # Plot 1: Boxplot of MP scores by hypoxia status
    mp_order <- c("MP1", "MP2", "MP3", "MP4", "MP5")
    mp_colors <- c("MP1" = "#0072B2", "MP2" = "#E69F00", "MP3" = "#009E73", "MP4" = "#CC79A7", "MP5" = "#D55E00")

    p1 <- ggplot(sample_data, aes(x = Metaprogram, y = Score, fill = hypoxia_status)) +
      geom_boxplot(alpha = 0.8, outlier.size = 1) +
      geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
      scale_fill_manual(
        values = c("Hypoxic" = "#E74C3C", "Normoxic" = "#3498DB"),
        name = "Region"
      ) +
      scale_x_discrete(limits = mp_order) +
      labs(
        x = "Metaprogram",
        y = "MP Score",
        title = "MP Score Distribution by Region"
      ) +
      theme_publication() +
      theme(legend.position = "right")
    
    plots[[1]] <- p1
    
    # Plot 2: Stacked bar chart of mean composition
    composition_long <- composition %>%
      pivot_longer(
        cols = -Metaprogram,
        names_to = "Region",
        values_to = "Mean_Score"
      )
    
    p2 <- ggplot(composition_long, aes(x = Region, y = Mean_Score, fill = Metaprogram)) +
      geom_col(alpha = 0.8, color = "black", size = 0.3) +
      scale_fill_manual(
        values = mp_colors,
        name = "Metaprogram"
      ) +
      labs(
        x = "Region Type",
        y = "Mean MP Score",
        title = "Mean MP Composition"
      ) +
      theme_publication() +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
    
    plots[[2]] <- p2
    
    # Plot 3: Violin plot
    p3 <- ggplot(sample_data, aes(x = Metaprogram, y = Score, fill = hypoxia_status)) +
      geom_violin(alpha = 0.8, position = position_dodge(width = 0.8)) +
      geom_point(
        aes(color = hypoxia_status),
        alpha = 0.2,
        size = 1,
        position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2)
      ) +
      scale_fill_manual(
        values = c("Hypoxic" = "#E74C3C", "Normoxic" = "#3498DB"),
        name = "Region"
      ) +
      scale_color_manual(
        values = c("Hypoxic" = "#C0392B", "Normoxic" = "#2980B9"),
        guide = "none"
      ) +
      scale_x_discrete(limits = mp_order) +
      labs(
        x = "Metaprogram",
        y = "MP Score",
        title = "MP Score Distribution (Violin)"
      ) +
      theme_publication()
    
    plots[[3]] <- p3
    
    # Combine plots
    combined_plot <- wrap_plots(plots, ncol = 3) +
      plot_annotation(
        title = paste("Sample:", sample),
        subtitle = "Metaprogram Distribution: Hypoxic vs Normoxic Regions",
        theme = theme(
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray50")
        )
      )
    
    # Save PDF
    pdf_filename <- paste0("Sample_", gsub("\\.", "_", sample), "_MP_Distribution.pdf")
    pdf_path <- file.path(OUTPUT_3B, pdf_filename)
    
    ggsave(pdf_path, combined_plot, width = 15, height = 5, dpi = 300)
    cat(" \n")
    
    generated_files_3b <- c(generated_files_3b, pdf_filename)
    
  }, error = function(e) {
    cat(" Error:", e$message, "\n")
  })
}

cat("\nTask 3B complete:", length(generated_files_3b), "sample plots generated\n\n")

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("TASK 3 COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Task 3A Summary:\n")
cat("  Generated:", length(generated_files_3a), "cell line comparative boxplots\n")
cat("  Output directory:", OUTPUT_3A, "\n")
if (length(generated_files_3a) > 0) {
  cat("  Files:\n")
  for (f in sort(generated_files_3a)) {
    cat("    ", f, "\n")
  }
}

cat("\nTask 3B Summary:\n")
cat("  Generated:", length(generated_files_3b), "individual sample MP distributions\n")
cat("  Output directory:", OUTPUT_3B, "\n")
if (length(generated_files_3b) > 5) {
  cat("  Files: (showing first 5)\n")
  for (f in sort(generated_files_3b)[1:5]) {
    cat("    ", f, "\n")
  }
  cat("    ... and", length(generated_files_3b) - 5, "more\n")
} else {
  cat("  Files:\n")
  for (f in sort(generated_files_3b)) {
    cat("    ", f, "\n")
  }
}

cat("\nAll Task 3 files ready for Task 4 reorganization\n")
