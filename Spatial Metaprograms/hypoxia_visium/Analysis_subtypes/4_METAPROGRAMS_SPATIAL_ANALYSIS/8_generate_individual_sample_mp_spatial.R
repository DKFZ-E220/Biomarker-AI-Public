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
# Task 2: Generate Individual Sample Metaprogram Spatial Visualizations
# ==============================================================================
# Purpose: Create spatial maps showing metaprogram distribution + hypoxia status
# for each sample individually, like the gene expression maps format.
# Output: One PDF per sample to designated output directory
# ==============================================================================

library(tidyverse)
library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Data paths
BASE_DIR <- file.path(dirname(dirname(analysis_subtypes_dir)))
MASTER_SPOT_TABLE <- file.path(analysis_subtypes_dir, "2_METAPROGRAMS_INTEGRATION", 
                                "data", "master_spot_table.csv")

# Output configuration
OUTPUT_PARENT_DIR <- file.path(BASE_DIR, "Writing")
PDF_OUTPUT_DIR <- file.path(OUTPUT_PARENT_DIR, "Task2_Individual_Sample_Spatial_MP")
dir.create(PDF_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Theme for publication quality
theme_spatial <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
      axis.title = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 9),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right",
      panel.grid = element_blank(),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
}

cat("═══════════════════════════════════════════════════════════════\n")
cat("   TASK 2: INDIVIDUAL SAMPLE METAPROGRAM SPATIAL MAPS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Output directory:", PDF_OUTPUT_DIR, "\n\n")

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading data...\n")

# Load master spot table with all data already computed
if (!file.exists(MASTER_SPOT_TABLE)) {
  stop("Master spot table not found:", MASTER_SPOT_TABLE)
}
master_data <- read_csv(MASTER_SPOT_TABLE, show_col_types = FALSE)
cat("Master data:", nrow(master_data), "spots\n")

# Get unique samples
all_samples <- sort(unique(master_data$sample))
cat("Samples:", length(all_samples), "\n")
cat("  ", paste(all_samples, collapse = ", "), "\n\n")

# ==============================================================================
# LOAD METAPROGRAMS
# ==============================================================================

cat("Setting up metaprograms...\n")
mp_names <- c("MP1", "MP2", "MP3", "MP4", "MP5")
cat("Metaprograms ready:", paste(mp_names, collapse = ", "), "\n")

# Check annotation column
has_annotation <- "pathologist_label" %in% colnames(master_data)
cat("Pathologist annotation available:", has_annotation, "\n\n")

# ==============================================================================
# CALCULATE GLOBAL MP SCALES (Tirosh approach)
# ==============================================================================

cat("Calculating global MP score ranges...\n")

# Calculate min/max for each MP across all samples for consistent scaling
mp_score_cols <- c("MP1", "MP2", "MP3", "MP4", "MP5")
mp_scales <- list()

for (mp in mp_score_cols) {
  mp_min <- min(master_data[[mp]], na.rm = TRUE)
  mp_max <- max(master_data[[mp]], na.rm = TRUE)
  # Guard against all-NA columns (returns Inf/-Inf) → use safe fallback
  if (!is.finite(mp_min) || !is.finite(mp_max) || mp_min == mp_max) {
    mp_min <- 0; mp_max <- 1
    cat("  ", mp, ": [fallback 0,1 — no finite scores]\n", sep = "")
  } else {
    cat("  ", mp, ": [", round(mp_min, 3), ", ", round(mp_max, 3), "]\n", sep = "")
  }
  mp_scales[[mp]] <- c(min = mp_min, max = mp_max)
}
cat("\n")

# ==============================================================================
# FUNCTION: Generate Individual Sample Spatial Plot
# ==============================================================================

generate_sample_spatial_plot <- function(sample_name, master_data, mp_scales) {
  
  # Filter data for this sample
  sample_data <- master_data %>% 
    filter(sample == !!sample_name)
  
  if (nrow(sample_data) == 0) {
    warning("No data found for", sample_name)
    return(NULL)
  }
  
  plot_data <- sample_data %>%
    mutate(
      hypoxia_status = case_when(
        hypoxia_label == 1 ~ "HYPOXIC",
        hypoxia_label == 0 ~ "NORMOXIC",
        TRUE ~ "Unknown"
      )
    )
  
  # Identify dominant metaprogram — guard against all-NA rows
  mp_score_cols <- c("MP1", "MP2", "MP3", "MP4", "MP5")
  mp_matrix <- as.matrix(plot_data[, mp_score_cols])
  dominant_mp_idx <- apply(mp_matrix, 1, function(x) {
    idx <- which.max(x)
    if (length(idx) == 0L) NA_integer_ else idx
  })
  plot_data$dominant_mp <- mp_score_cols[dominant_mp_idx]
  
  # Create plots
  plots <- list()
  
  # MP score panels (Wong palette, global scale)
  mp_info <- list(
    MP1 = list(color = "#0072B2", label = "MP1 — EMT/Mesenchymal"),
    MP2 = list(color = "#E69F00", label = "MP2 — Hypoxia Response"),
    MP3 = list(color = "#009E73", label = "MP3 — Squamous Diff."),
    MP4 = list(color = "#CC79A7", label = "MP4 — Fibroblast/Stromal"),
    MP5 = list(color = "#D55E00", label = "MP5 — Keratinocyte/Basal")
  )
  for (i in seq_along(mp_info)) {
    mp  <- names(mp_info)[i]
    col <- mp_info[[mp]]$color
    lbl <- mp_info[[mp]]$label
    p <- ggplot(plot_data, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
                               color = .data[[mp]])) +
      geom_point(size = 2.5, alpha = 0.8, na.rm = TRUE) +
      scale_color_gradient(low = "white", high = col, name = "Score",
                           limits = mp_scales[[mp]], oob = scales::squish) +
      scale_y_reverse() +
      coord_equal() +
      labs(title = lbl, x = "X", y = "Y") +
      theme_spatial()
    plots[[i]] <- p
  }

  # Plot 6: Hypoxia Status
  plots[[6]] <- ggplot(plot_data, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
                                      color = hypoxia_status)) +
    geom_point(size = 2.5, alpha = 0.8, na.rm = TRUE) +
    scale_color_manual(values = c("HYPOXIC" = "#D55E00", "NORMOXIC" = "#0072B2",
                                  "Unknown" = "gray80"), name = "Status",
                       na.value = "gray80") +
    scale_y_reverse() +
    coord_equal() +
    labs(title = "Hypoxia Status (PIMO)", x = "X", y = "Y") +
    theme_spatial()

  # Plot 7: Dominant Metaprogram
  plots[[7]] <- ggplot(plot_data, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
                                      color = dominant_mp)) +
    geom_point(size = 2.5, alpha = 0.8, na.rm = TRUE) +
    scale_color_manual(
      values = c("MP1"="#0072B2","MP2"="#E69F00","MP3"="#009E73","MP4"="#CC79A7","MP5"="#D55E00"),
      name = "Dominant MP", na.value = "gray80") +
    scale_y_reverse() +
    coord_equal() +
    labs(title = "Dominant Metaprogram", x = "X", y = "Y") +
    theme_spatial()

  # Plot 8: Pathologist annotation (if available)
  if (has_annotation && "pathologist_label" %in% colnames(plot_data)) {
    annot_colors <- c(
      "Tumor"                   = "#E69F00",
      "Necrosis"                = "#000000",
      "Stroma"                  = "#56B4E9",
      "Keratinization"          = "#009E73",
      "Blood vessels"           = "#D55E00",
      "Inflammatory infiltrate" = "#CC79A7",
      "Immune cells"            = "#CC79A7",
      "Edema"                   = "#0072B2",
      "Muscle"                  = "#8B4513",
      "Epidermis"               = "#F0E442",
      "Bacteria"                = "#999999",
      "Ambiguous"               = "#AAAAAA",  # border spots: darker grey, distinct from Unannotated
      "Unannotated"             = "#DDDDDD"   # no annotation drawn: light grey
    )
    plots[[8]] <- ggplot(plot_data, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres,
                                        color = pathologist_label)) +
      geom_point(size = 2.5, alpha = 0.8, na.rm = TRUE) +
      scale_color_manual(values = annot_colors, na.value = "#DDDDDD",
                         name = "Pathologist") +
      scale_y_reverse() +
      coord_equal() +
      labs(title = "Pathologist Annotation", x = "X", y = "Y") +
      theme_spatial()
  }

  return(list(plots = plots, data = plot_data, sample = sample_name))
}

# ==============================================================================
# GENERATE PLOTS FOR ALL SAMPLES
# ==============================================================================

cat("Generating spatial plots...\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

generated_files <- character()
success_samples <- character()
failed_samples <- character()

for (sample in all_samples) {
  
  tryCatch({
    cat("Processing", sample, "...")
    
    result <- generate_sample_spatial_plot(sample, master_data, mp_scales)
    
    if (is.null(result)) {
      failed_samples <- c(failed_samples, sample)
      cat(" Failed\n")
      next
    }
    
    # 8 panels: MP1-5 + Hypoxia + Dominant + Annotation (2 rows × 4 cols)
    combined_plot <- wrap_plots(result$plots, ncol = 4, nrow = 2) +
      plot_annotation(
        title = paste("Sample:", sample),
        subtitle = "Metaprogram Spatial Distribution + Hypoxia Status",
        theme = theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray50")
        )
      )
    
    # Save PDF
    pdf_filename <- paste0("Spatial_MP_", gsub("\\.", "_", sample), ".pdf")
    pdf_path <- file.path(PDF_OUTPUT_DIR, pdf_filename)
    
    ggsave(pdf_path, combined_plot, width = 20, height = 10, dpi = 300)
    cat(" Saved\n")
    
    generated_files <- c(generated_files, pdf_filename)
    success_samples <- c(success_samples, sample)
    
  }, error = function(e) {
    cat(" Error:", e$message, "\n")
    failed_samples <<- c(failed_samples, sample)
  })
}

cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("TASK 2 COMPLETE\n\n")

cat("Summary:\n")
cat("  Total samples processed:", length(all_samples), "\n")
cat("  Successful:", length(success_samples), "\n")
cat("  Failed:", length(failed_samples), "\n")
cat("  PDFs generated:", length(generated_files), "\n\n")

if (length(generated_files) > 0) {
  cat("Generated files:\n")
  for (f in sort(generated_files)) {
    cat("  ", f, "\n")
  }
}

if (length(failed_samples) > 0) {
  cat("\nFailed samples:\n")
  for (s in failed_samples) {
    cat("  ", s, "\n")
  }
}

cat("\nOutput directory:", PDF_OUTPUT_DIR, "\n\n")

cat("Each PDF contains 6 spatial maps:\n")
cat("  1. MP1 (Hypoxia Response) - viridis\n")
cat("  2. MP2 (Squamous Differentiation) - plasma\n")
cat("  3. MP3 (Baseline Epithelial) - inferno\n")
cat("  4. MP4 (EMT-Invasive) - magma\n")
cat("  5. Hypoxia Status (Red=Hypoxic, Blue=Normoxic)\n")
cat("  6. Dominant Metaprogram classification\n\n")

cat("Task 2 files ready for Task 4 reorganization\n")
