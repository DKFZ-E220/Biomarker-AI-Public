#!/usr/bin/env Rscript

# Resolve script directory
script_dir <- tryCatch({
  frames <- sys.frames(); ofile <- NULL
  for (f in rev(frames)) { if (!is.null(f$ofile)) { ofile <- f$ofile; break } }
  if (!is.null(ofile) && nchar(ofile) > 0) dirname(normalizePath(ofile)) else getwd()
}, error = function(e) getwd())
setwd(script_dir)

analysis_subtypes_dir <- dirname(script_dir)
spatial_data_dir      <- dirname(dirname(analysis_subtypes_dir))

# ==============================================================================
# SPOT-LEVEL SUBTYPE vs METAPROGRAM — Side-by-Side Spatial Comparison
# ==============================================================================
# For each sample:
#   Left panel  — spot assigned to bulk molecular subtype (centroid correlation)
#   Right panel — spot assigned to metaprogram (NMF scores)
# Output: one PDF per sample + combined multi-page PDF
# ==============================================================================

library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(readxl)

# ==============================================================================
# PATHS
# ==============================================================================

QC_DIR            <- file.path(analysis_subtypes_dir, "0_QC_AND_PREPROCESSING", "data")
CENTROID_FILE     <- file.path(spatial_data_dir, "Bulk data", "updated_centroid_values.xlsx")
MASTER_SPOT_TABLE <- file.path(analysis_subtypes_dir, "2_METAPROGRAMS_INTEGRATION",
                                "data", "master_spot_table.csv")
OUTPUT_DIR        <- file.path(spatial_data_dir, "Writing", "SpotSubtype_vs_Metaprogram")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# COLOUR PALETTES
# ==============================================================================

subtype_colors <- c(
  Atypical     = "#00ADB5",
  Basal        = "#56B4E9",
  Mesenchymal  = "#D55E00",
  Classical    = "#F9ED69",
  Heterogenous = "#CCCCCC",  # grey — unresolved spots
  Unknown      = "#FFFFFF"
)

mp_colors <- c(
  MP1 = "#0072B2",   # EMT/Mesenchymal
  MP2 = "#E69F00",   # Hypoxia Response
  MP3 = "#009E73",   # Squamous Differentiation
  MP4 = "#CC79A7",   # Fibroblast/Stromal
  MP5 = "#D55E00"    # Keratinocyte/Basal
)

mp_labels <- c(
  MP1 = "MP1 — EMT/Mesenchymal",
  MP2 = "MP2 — Hypoxia Response",
  MP3 = "MP3 — Squamous Diff.",
  MP4 = "MP4 — Fibroblast/Stromal",
  MP5 = "MP5 — Keratinocyte/Basal"
)

# ==============================================================================
# LOAD SHARED DATA
# ==============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("   SPOT-LEVEL SUBTYPE vs METAPROGRAM COMPARISON\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Loading centroid values...\n")
centroid_values <- read_excel(CENTROID_FILE, sheet = "Sheet1") %>%
  as.data.frame()
rownames(centroid_values) <- centroid_values$Gene_Symbol
centroid_values <- centroid_values[, -1]
cat("Centroids:", ncol(centroid_values), "subtypes,", nrow(centroid_values), "genes\n")
cat("  Subtypes:", paste(colnames(centroid_values), collapse = ", "), "\n\n")

cat("Loading master spot table (MP assignments)...\n")
master_spots <- read_csv(MASTER_SPOT_TABLE, show_col_types = FALSE)
cat("Spots:", nrow(master_spots), "\n\n")

# Samples with MP data
samples_with_mp <- sort(unique(master_spots$sample))
cat("Samples to process:", length(samples_with_mp), "\n")
cat(" ", paste(samples_with_mp, collapse = "\n  "), "\n\n")

# ==============================================================================
# MEDIAN CENTERING FUNCTION
# ==============================================================================

median_center <- function(x) x - median(x, na.rm = TRUE)

# ==============================================================================
# PROCESS EACH SAMPLE
# ==============================================================================

cat("Processing samples...\n")
cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")

pdf_pages        <- list()
failed           <- character()
spot_subtype_all <- list()   # collector for pooled summary

for (sample_name in samples_with_mp) {

  cat("▶ Sample:", sample_name, "\n")

  tryCatch({

    # ── Load QC Seurat object ────────────────────────────────────────────────
    rds_file <- file.path(QC_DIR, paste0(sample_name, "_qc_with_metadata.rds"))
    if (!file.exists(rds_file)) {
      cat("  RDS not found:", rds_file, "\n\n")
      failed <- c(failed, sample_name)
      next
    }
    seu <- readRDS(rds_file)
    cat("  Loaded:", ncol(seu), "spots\n")

    # ── Normalize if needed, then median-center ──────────────────────────────
    seu <- NormalizeData(seu, assay = "Spatial",
                         normalization.method = "LogNormalize", verbose = FALSE)
    expr_mat <- as.matrix(GetAssayData(seu, assay = "Spatial", layer = "data"))
    expr_mat  <- t(apply(expr_mat, 1, median_center))

    # ── Centroid correlation ─────────────────────────────────────────────────
    common_genes <- intersect(rownames(expr_mat), rownames(centroid_values))
    cat("  Common genes for centroid correlation:", length(common_genes), "\n")

    expr_sub     <- expr_mat[common_genes, , drop = FALSE]
    centroid_sub <- as.matrix(centroid_values[common_genes, , drop = FALSE])

    # Remove all-zero rows
    keep <- rowSums(expr_sub != 0) > 0
    expr_sub     <- expr_sub[keep, , drop = FALSE]
    centroid_sub <- centroid_sub[keep, , drop = FALSE]

    # Correlate each spot vs each centroid
    cor_mat <- cor(expr_sub, centroid_sub, method = "pearson")  # spots × subtypes

    # Assign best subtype per spot with post-processing (same logic as bulk)
    # A spot gets a definitive subtype only if:
    #   max_cor > 0.2  AND  (1st - 2nd) > 0.2  → otherwise "Heterogenous"
    max_cor <- apply(cor_mat, 1, function(r) max(r, na.rm = TRUE))

    assigned_idx <- apply(cor_mat, 1, function(r) {
      sorted <- sort(r, decreasing = TRUE)
      if (!is.na(sorted[1]) && sorted[1] > 0.2 && (sorted[1] - sorted[2]) > 0.2) {
        which.max(r)
      } else {
        NA_integer_
      }
    })

    best_subtype <- ifelse(is.na(assigned_idx),
                           "Heterogenous",
                           colnames(centroid_sub)[assigned_idx])

    seu@meta.data$BulkSubtype    <- best_subtype
    seu@meta.data$MaxCorrelation <- max_cor
    cat("  Subtype distribution:",
        paste(names(table(best_subtype)), table(best_subtype), sep = "=", collapse = ", "), "\n")

    # ── Add MP assignments from master spot table ────────────────────────────
    mp_df <- master_spots %>%
      filter(sample == sample_name) %>%
      select(barcode, dominant_metaprogram)

    # Robust join — handles barcode mismatches between QC object and master table
    meta_mp <- data.frame(barcode = colnames(seu), stringsAsFactors = FALSE) %>%
      left_join(mp_df, by = "barcode") %>%
      mutate(MetaProgram = ifelse(is.na(dominant_metaprogram), "Unassigned",
                                  dominant_metaprogram))

    seu@meta.data$MetaProgram <- meta_mp$MetaProgram
    cat("  MP distribution:",
        paste(names(table(seu@meta.data$MetaProgram)),
              table(seu@meta.data$MetaProgram), sep = "=", collapse = ", "), "\n")

    # ── Side-by-side spatial plots ───────────────────────────────────────────
    # Left: bulk subtype per spot
    p_subtype <- SpatialDimPlot(
      seu,
      group.by = "BulkSubtype",
      cols     = subtype_colors,
      pt.size.factor = 1.3,
      image.alpha    = 0.6
    ) +
      labs(title    = "Bulk Subtype (Spot-Level)",
           subtitle = sample_name,
           fill     = "Subtype") +
      theme(plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 9,  hjust = 0.5, colour = "grey40"),
            legend.title  = element_text(size = 9,  face = "bold"),
            legend.text   = element_text(size = 8))

    # Right: metaprogram per spot
    mp_col_use <- mp_colors[names(mp_colors) %in% unique(seu$MetaProgram)]
    p_mp <- SpatialDimPlot(
      seu,
      group.by = "MetaProgram",
      cols     = mp_col_use,
      pt.size.factor = 1.3,
      image.alpha    = 0.6
    ) +
      labs(title    = "Metaprogram (NMF)",
           subtitle = sample_name,
           fill     = "Metaprogram") +
      theme(plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 9,  hjust = 0.5, colour = "grey40"),
            legend.title  = element_text(size = 9,  face = "bold"),
            legend.text   = element_text(size = 8))

    combined <- p_subtype + p_mp +
      plot_annotation(
        title = sample_name,
        theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
      )

    # ── Collect spot-level subtype counts for pooled summary ────────────────
    spot_subtype_all[[sample_name]] <- data.frame(
      sample   = sample_name,
      subtype  = best_subtype,
      stringsAsFactors = FALSE
    )

    # Save individual PDF
    out_file <- file.path(OUTPUT_DIR, paste0(sample_name, "_subtype_vs_mp.pdf"))
    ggsave(out_file, combined, width = 14, height = 7, dpi = 300)
    cat("  Saved:", basename(out_file), "\n\n")

    pdf_pages[[sample_name]] <- combined

    # Free memory
    rm(seu, expr_mat, expr_sub, cor_mat)
    gc(verbose = FALSE)

  }, error = function(e) {
    cat("  Error:", e$message, "\n\n")
    failed <<- c(failed, sample_name)
  })
}

# ==============================================================================
# COMBINED MULTI-PAGE PDF
# ==============================================================================

cat("Writing combined PDF...\n")
combined_pdf <- file.path(OUTPUT_DIR, "ALL_SAMPLES_subtype_vs_mp.pdf")
pdf(combined_pdf, width = 14, height = 7)
for (nm in names(pdf_pages)) {
  tryCatch(print(pdf_pages[[nm]]), error = function(e) cat("  Skip", nm, "\n"))
}
dev.off()
cat("  Saved: ALL_SAMPLES_subtype_vs_mp.pdf\n\n")

# ==============================================================================
# POOLED SPOT-LEVEL SUBTYPE SUMMARY — CSV + BAR FIGURE
# ==============================================================================

cat("Writing pooled spot-level subtype summary...\n")

if (length(spot_subtype_all) > 0) {

  pooled <- bind_rows(spot_subtype_all)

  # Per-sample counts
  per_sample <- pooled %>%
    count(sample, subtype) %>%
    group_by(sample) %>%
    mutate(total = sum(n),
           pct   = round(n / total * 100, 2)) %>%
    ungroup()

  # Overall pooled counts across all spots
  overall <- pooled %>%
    count(subtype) %>%
    mutate(total = sum(n),
           pct   = round(n / total * 100, 2)) %>%
    arrange(desc(n))

  # Write CSVs
  write.csv(per_sample,
            file.path(OUTPUT_DIR, "Pooled_SpotSubtype_PerSample.csv"),
            row.names = FALSE)
  write.csv(overall,
            file.path(OUTPUT_DIR, "Pooled_SpotSubtype_Overall.csv"),
            row.names = FALSE)
  cat("  Saved: Pooled_SpotSubtype_PerSample.csv\n")
  cat("  Saved: Pooled_SpotSubtype_Overall.csv\n")
  cat("  Overall distribution:\n")
  for (i in seq_len(nrow(overall))) {
    cat(sprintf("    %-14s  %6d spots  (%5.1f%%)\n",
                overall$subtype[i], overall$n[i], overall$pct[i]))
  }
  cat("\n")

  # Bar figure — overall % per subtype
  bar_colors_use <- subtype_colors[names(subtype_colors) %in% overall$subtype]

  fig_bar <- ggplot(overall,
                    aes(x = reorder(subtype, -pct), y = pct, fill = subtype)) +
    geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%", pct)),
              vjust = -0.4, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = bar_colors_use) +
    scale_y_continuous(limits = c(0, max(overall$pct) * 1.15),
                       labels = function(x) paste0(x, "%")) +
    labs(title    = "Spot-Level Bulk Subtype Assignment (All Samples Pooled)",
         subtitle = sprintf("%d spots across %d samples",
                            unique(overall$total), length(spot_subtype_all)),
         x = "Assigned Bulk Subtype",
         y = "% of All Spots") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, colour = "grey40"),
      axis.title    = element_text(size = 11, face = "bold"),
      axis.text     = element_text(size = 10),
      legend.position = "none",
      panel.grid.minor  = element_blank(),
      panel.grid.major.x = element_blank()
    )

  ggsave(file.path(OUTPUT_DIR, "Fig_SpotSubtype_Overall_Bar.pdf"),
         fig_bar, width = 7, height = 5, dpi = 300)
  cat("  Saved: Fig_SpotSubtype_Overall_Bar.pdf\n\n")

} else {
  cat("  No spot data collected — skipping pooled summary\n\n")
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("DONE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
cat("Processed:", length(pdf_pages), "samples\n")
if (length(failed) > 0) cat("Failed:   ", paste(failed, collapse = ", "), "\n")
cat("Output:", OUTPUT_DIR, "\n\n")
