#!/usr/bin/env Rscript

# ==============================================================================
# MANUAL QC INSPECTION — Per-Sample Spatial QC Plots
# ==============================================================================
# For each sample generates a PDF with 4 spatial panels:
#   1. Total UMI counts       — low = off-tissue / mouse / necrotic
#   2. Genes detected         — low = poor quality / mouse contamination
#   3. Mitochondrial %        — high = necrotic / dying cells
#   4. UMI histogram          — bimodal = two populations (real tissue vs noise)
#
# Review PDFs manually → fill in per_sample_thresholds.csv → run apply_qc_filters.R
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# ==============================================================================
# PATHS
# ==============================================================================

source(here::here("config.R"))  # defines DATA_DIR, RESULTS_DIR, ANALYSIS_SUBTYPES_DIR
script_dir    <- file.path(ANALYSIS_SUBTYPES_DIR, "0_QC_AND_PREPROCESSING")
qc_data_dir   <- file.path(script_dir, "data")
output_dir    <- file.path(script_dir, "QC_inspection")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

rds_files <- list.files(qc_data_dir,
                         pattern = "_qc_with_metadata\\.rds$",
                         full.names = TRUE)

cat("Found", length(rds_files), "samples\n\n")

# ==============================================================================
# PROCESS EACH SAMPLE
# ==============================================================================

summary_rows <- list()

for (rds_file in rds_files) {

  sample_name <- gsub("_qc_with_metadata\\.rds$", "", basename(rds_file))
  cat("Processing:", sample_name, "\n")

  seu <- tryCatch(readRDS(rds_file), error = function(e) {
    cat("  ERROR loading:", e$message, "\n"); NULL
  })
  if (is.null(seu)) next

  # ── Mitochondrial % ─────────────────────────────────────────────────────────
  mt_genes <- grep("^MT-", rownames(seu), value = TRUE, ignore.case = TRUE)
  if (length(mt_genes) > 0) {
    seu[["pct_mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  } else {
    seu[["pct_mt"]] <- 0
    cat("  No MT genes found — pct_mt set to 0\n")
  }

  # ── Coordinates ─────────────────────────────────────────────────────────────
  coords_df <- GetTissueCoordinates(seu)
  meta      <- seu@meta.data
  meta$barcode <- rownames(meta)
  meta$x       <- coords_df$x[match(rownames(meta), rownames(coords_df))]
  meta$y       <- coords_df$y[match(rownames(meta), rownames(coords_df))]

  # ── Summary stats ────────────────────────────────────────────────────────────
  summary_rows[[sample_name]] <- data.frame(
    sample       = sample_name,
    n_spots      = ncol(seu),
    median_UMI   = round(median(meta$nCount_Spatial)),
    median_genes = round(median(meta$nFeature_Spatial)),
    median_pct_mt= round(median(meta$pct_mt), 1),
    pct_low_UMI  = round(mean(meta$nCount_Spatial < 500) * 100, 1)
  )

  # ── PANEL 1: UMI counts ──────────────────────────────────────────────────────
  p1 <- ggplot(meta, aes(x = x, y = y, colour = log1p(nCount_Spatial))) +
    geom_point(size = 0.6, alpha = 0.8) +
    scale_colour_viridis_c(name = "log1p(UMI)", option = "C") +
    coord_equal() +
    labs(title = "Total UMI counts",
         subtitle = paste0("median = ", median(meta$nCount_Spatial))) +
    theme_void() +
    theme(plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
          plot.subtitle = element_text(size = 9, hjust = 0.5, colour = "grey40"),
          legend.title  = element_text(size = 8),
          legend.text   = element_text(size = 7))

  # ── PANEL 2: Genes detected ──────────────────────────────────────────────────
  p2 <- ggplot(meta, aes(x = x, y = y, colour = log1p(nFeature_Spatial))) +
    geom_point(size = 0.6, alpha = 0.8) +
    scale_colour_viridis_c(name = "log1p(genes)", option = "D") +
    coord_equal() +
    labs(title = "Genes detected",
         subtitle = paste0("median = ", median(meta$nFeature_Spatial))) +
    theme_void() +
    theme(plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
          plot.subtitle = element_text(size = 9, hjust = 0.5, colour = "grey40"),
          legend.title  = element_text(size = 8),
          legend.text   = element_text(size = 7))

  # ── PANEL 3: Mitochondrial % ─────────────────────────────────────────────────
  p3 <- ggplot(meta, aes(x = x, y = y, colour = pct_mt)) +
    geom_point(size = 0.6, alpha = 0.8) +
    scale_colour_gradient(low = "grey90", high = "red",
                          name = "% MT") +
    coord_equal() +
    labs(title = "Mitochondrial %",
         subtitle = paste0("median = ", round(median(meta$pct_mt), 1), "%  |  ",
                           "high MT (>25%) = ",
                           sum(meta$pct_mt > 25), " spots")) +
    theme_void() +
    theme(plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
          plot.subtitle = element_text(size = 9, hjust = 0.5, colour = "grey40"),
          legend.title  = element_text(size = 8),
          legend.text   = element_text(size = 7))

  # ── PANEL 4: UMI histogram ───────────────────────────────────────────────────
  p4 <- ggplot(meta, aes(x = log1p(nCount_Spatial))) +
    geom_histogram(bins = 80, fill = "#0072B2", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = log1p(500),  colour = "red",    linetype = "dashed") +
    geom_vline(xintercept = log1p(1000), colour = "orange", linetype = "dashed") +
    geom_vline(xintercept = log1p(2000), colour = "green",  linetype = "dashed") +
    annotate("text", x = log1p(500),  y = Inf, label = "500",
             vjust = 2, hjust = -0.1, size = 3, colour = "red") +
    annotate("text", x = log1p(1000), y = Inf, label = "1000",
             vjust = 2, hjust = -0.1, size = 3, colour = "orange") +
    annotate("text", x = log1p(2000), y = Inf, label = "2000",
             vjust = 2, hjust = -0.1, size = 3, colour = "darkgreen") +
    labs(title = "UMI distribution",
         x = "log1p(UMI)", y = "Spots") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

  # ── Combine & save ───────────────────────────────────────────────────────────
  combined <- (p1 + p2) / (p3 + p4) +
    plot_annotation(
      title    = sample_name,
      subtitle = paste0("Total spots: ", ncol(seu),
                        "  |  Median UMI: ", median(meta$nCount_Spatial),
                        "  |  Median genes: ", median(meta$nFeature_Spatial),
                        "  |  Median MT%: ", round(median(meta$pct_mt), 1)),
      theme = theme(
        plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, colour = "grey40")
      )
    )

  out_pdf <- file.path(output_dir, paste0(sample_name, "_QC_inspection.pdf"))
  ggsave(out_pdf, combined, width = 14, height = 10, dpi = 200)
  cat("  Saved:", basename(out_pdf), "\n\n")

  rm(seu); gc(verbose = FALSE)
}

# ==============================================================================
# SAVE SUMMARY TABLE + THRESHOLD TEMPLATE
# ==============================================================================

summary_df <- bind_rows(summary_rows)
write.csv(summary_df,
          file.path(output_dir, "QC_summary_all_samples.csv"),
          row.names = FALSE)

# Template CSV for user to fill in per-sample thresholds
threshold_template <- summary_df %>%
  select(sample) %>%
  mutate(
    min_UMI    = 500,   # adjust per sample after reviewing histograms
    max_pct_mt = 25,    # spots above this % MT removed (necrotic)
    notes      = ""     # e.g. "bimodal at 1500" or "necrotic patch top-right"
  )

write.csv(threshold_template,
          file.path(output_dir, "per_sample_thresholds.csv"),
          row.names = FALSE)

cat("═══════════════════════════════════════════════════════\n")
cat("DONE\n")
cat("═══════════════════════════════════════════════════════\n\n")
cat("Review PDFs in:  ", output_dir, "\n")
cat("Fill thresholds: per_sample_thresholds.csv\n")
cat("Then run:        apply_qc_filters.R\n\n")
cat("Summary table:\n")
print(summary_df)
