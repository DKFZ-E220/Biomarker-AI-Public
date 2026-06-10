#!/usr/bin/env Rscript

# Resolve script directory
script_dir <- tryCatch({
  frames <- sys.frames()
  ofile <- NULL
  for (f in rev(frames)) { if (!is.null(f$ofile)) { ofile <- f$ofile; break } }
  if (!is.null(ofile) && nchar(ofile) > 0) dirname(normalizePath(ofile)) else getwd()
}, error = function(e) getwd())
setwd(script_dir)
analysis_subtypes_dir <- if (grepl("4_METAPROGRAMS_SPATIAL_ANALYSIS", getwd(), fixed = TRUE)) dirname(getwd()) else getwd()

# ==============================================================================
# ANNOTATION × METAPROGRAM × HYPOXIA COMPARATIVE PLOTS
# ==============================================================================
# Connects pathologist H&E annotations with NMF metaprograms and hypoxia status.
# Requires pathologist_label column in master_spot_table (run script 1 first).
#
# Outputs (Writing/Annotation_MP_Hypoxia/):
#   1. MP_Score_by_Tissue_Type.pdf       — boxplots per MP across tissue labels
#   2. Hypoxia_by_Tissue_Type.pdf        — % hypoxic spots per tissue type
#   3. Dominant_MP_per_Tissue.pdf        — stacked bars: tissue → MP composition
#   4. Tissue_per_DominantMP.pdf         — stacked bars: dominant MP → tissue composition
#   5. MP_Heatmap_Tissue_Hypoxia.pdf     — heatmap: mean MP score × tissue × hypoxia
#   6. Tumor_MP_Hypoxia_Shift.pdf        — within-tumor MP shift hypoxic vs normoxic
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
})

# ==============================================================================
# PATHS
# ==============================================================================

BASE_DIR         <- dirname(dirname(analysis_subtypes_dir))
MASTER_SPOT_TABLE <- file.path(analysis_subtypes_dir, "2_METAPROGRAMS_INTEGRATION",
                                "data", "master_spot_table.csv")
OUTPUT_DIR <- file.path(BASE_DIR, "Writing", "Annotation_MP_Hypoxia")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("═══════════════════════════════════════════════════════════════\n")
cat("   ANNOTATION × METAPROGRAM × HYPOXIA COMPARATIVE PLOTS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
cat("Output:", OUTPUT_DIR, "\n\n")

# ==============================================================================
# LOAD DATA
# ==============================================================================

if (!file.exists(MASTER_SPOT_TABLE)) stop("master_spot_table.csv not found: ", MASTER_SPOT_TABLE)
master <- read_csv(MASTER_SPOT_TABLE, show_col_types = FALSE)
cat("Loaded:", nrow(master), "spots\n")

if (!"pathologist_label" %in% colnames(master)) {
  stop("pathologist_label column not found — re-run script 1 (1_MP_visualization_hypoxia_CLEAN.Rmd) first.")
}

# ==============================================================================
# PALETTE & LABELS
# ==============================================================================

mp_colors <- c(
  MP1 = "#0072B2",  # EMT/Mesenchymal
  MP2 = "#E69F00",  # Hypoxia Response
  MP3 = "#009E73",  # Squamous Differentiation
  MP4 = "#CC79A7",  # Fibroblast/Stromal
  MP5 = "#D55E00"   # Keratinocyte/Basal
)
mp_labels <- c(
  MP1 = "MP1\nEMT/Mesen.",
  MP2 = "MP2\nHypoxia",
  MP3 = "MP3\nSquamous",
  MP4 = "MP4\nFibroblast",
  MP5 = "MP5\nKeratinocyte"
)
hypoxia_colors <- c("Hypoxic" = "#D55E00", "Normoxic" = "#0072B2")

annot_colors <- c(
  "Tumor"                   = "#E69F00",
  "Necrosis"                = "#111111",
  "Stroma"                  = "#56B4E9",
  "Keratinization"          = "#009E73",
  "Blood vessels"           = "#D55E00",
  "Inflammatory infiltrate" = "#CC79A7",
  "Immune cells"            = "#AA77A7",
  "Edema"                   = "#0072B2",
  "Muscle"                  = "#8B4513",
  "Epidermis"               = "#F0E442",
  "Bacteria"                = "#999999",
  "Ambiguous"               = "#AAAAAA",  # border spots: darker grey
  "Unannotated"             = "#DDDDDD"   # no annotation drawn: light grey
)

# Priority order for tissue types on axes (most important first)
tissue_order <- c("Tumor", "Necrosis", "Keratinization", "Stroma",
                  "Blood vessels", "Inflammatory infiltrate", "Immune cells",
                  "Edema", "Epidermis", "Muscle", "Bacteria", "Ambiguous", "Unannotated")

# ==============================================================================
# PREPARE DATA — filter out NA rows
# ==============================================================================

df <- master %>%
  mutate(
    hypoxia_status = case_when(
      hypoxia_label == 1 ~ "Hypoxic",
      hypoxia_label == 0 ~ "Normoxic",
      TRUE               ~ NA_character_
    ),
    pathologist_label = factor(pathologist_label, levels = tissue_order)
  ) %>%
  filter(!is.na(pathologist_label),
         !pathologist_label %in% c("Unannotated", "Ambiguous"),
         !is.na(dominant_metaprogram))

cat("After filtering (no Unannotated/Ambiguous/NA):", nrow(df), "spots\n")
cat("  Tissue types present:", paste(sort(unique(as.character(df$pathologist_label))), collapse=", "), "\n\n")

# MP score long format
mp_long <- df %>%
  pivot_longer(cols = c(MP1, MP2, MP3, MP4, MP5),
               names_to = "mp", values_to = "score") %>%
  filter(!is.na(score))

# Publication theme
theme_pub <- function(angle_x = 45) {
  theme_minimal() +
    theme(
      plot.title      = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle   = element_text(size = 10, hjust = 0.5, color = "gray40"),
      axis.title      = element_text(size = 11, face = "bold"),
      axis.text       = element_text(size = 9),
      axis.text.x     = element_text(angle = angle_x, hjust = 1),
      legend.title    = element_text(size = 10, face = "bold"),
      legend.text     = element_text(size = 9),
      panel.grid.major = element_line(color = "gray93", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      strip.text      = element_text(size = 10, face = "bold")
    )
}

# ==============================================================================
# FIGURE 1: MP Score Distributions by Tissue Type
# ==============================================================================

cat("Figure 1: MP scores by tissue type...\n")

p1 <- ggplot(mp_long, aes(x = pathologist_label, y = score, fill = mp)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3,
               position = position_dodge(width = 0.8), linewidth = 0.3) +
  scale_fill_manual(values = mp_colors, name = "Metaprogram",
                    labels = c(MP1="MP1 EMT", MP2="MP2 Hypoxia",
                               MP3="MP3 Squamous", MP4="MP4 Fibroblast", MP5="MP5 Keratinocyte")) +
  facet_wrap(~ mp, ncol = 5, labeller = labeller(mp = mp_labels)) +
  labs(title = "Metaprogram Score Distribution by Tissue Type",
       subtitle = "H&E pathologist annotation vs NMF metaprogram scores",
       x = "Tissue Type", y = "MP Score") +
  theme_pub() +
  theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "MP_Score_by_Tissue_Type.pdf"),
       p1, width = 18, height = 6, dpi = 300)
cat("  MP_Score_by_Tissue_Type.pdf\n")

# ==============================================================================
# FIGURE 2: % Hypoxic Spots per Tissue Type
# ==============================================================================

cat("Figure 2: Hypoxia fraction per tissue type...\n")

hypoxia_by_tissue <- df %>%
  filter(!is.na(hypoxia_status)) %>%
  group_by(pathologist_label) %>%
  summarise(
    n_total   = n(),
    n_hypoxic = sum(hypoxia_status == "Hypoxic"),
    pct_hypoxic = 100 * n_hypoxic / n_total,
    .groups = "drop"
  )

p2 <- ggplot(hypoxia_by_tissue, aes(x = pathologist_label, y = pct_hypoxic,
                                     fill = pathologist_label)) +
  geom_col(color = "black", linewidth = 0.3) +
  geom_text(aes(label = paste0(round(pct_hypoxic, 1), "%\n(n=", n_total, ")")),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = annot_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Hypoxia (PIMO) Frequency by Tissue Type",
       subtitle = "% spots classified as hypoxic within each pathologist annotation",
       x = "Tissue Type", y = "% Hypoxic Spots") +
  theme_pub()

ggsave(file.path(OUTPUT_DIR, "Hypoxia_by_Tissue_Type.pdf"),
       p2, width = 10, height = 6, dpi = 300)
cat("  Hypoxia_by_Tissue_Type.pdf\n")

# ==============================================================================
# FIGURE 3: Dominant MP Composition Within Each Tissue Type
# ==============================================================================

cat("Figure 3: MP composition per tissue type...\n")

mp_per_tissue <- df %>%
  group_by(pathologist_label, dominant_metaprogram) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(pathologist_label) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

p3 <- ggplot(mp_per_tissue, aes(x = pathologist_label, y = pct,
                                  fill = dominant_metaprogram)) +
  geom_col(color = "black", linewidth = 0.3) +
  scale_fill_manual(values = mp_colors, name = "Dominant MP",
                    labels = c(MP1="MP1 EMT/Mesen.", MP2="MP2 Hypoxia",
                               MP3="MP3 Squamous", MP4="MP4 Fibroblast",
                               MP5="MP5 Keratinocyte")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(title = "Metaprogram Composition Within Each Tissue Type",
       subtitle = "Dominant metaprogram assignment per pathologist label",
       x = "Tissue Type", y = "% Spots") +
  theme_pub()

ggsave(file.path(OUTPUT_DIR, "Dominant_MP_per_Tissue.pdf"),
       p3, width = 12, height = 6, dpi = 300)
cat("  Dominant_MP_per_Tissue.pdf\n")

# ==============================================================================
# FIGURE 4: Tissue Composition Within Each Dominant MP
# ==============================================================================

cat("Figure 4: Tissue composition per dominant MP...\n")

tissue_per_mp <- df %>%
  group_by(dominant_metaprogram, pathologist_label) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(dominant_metaprogram) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

p4 <- ggplot(tissue_per_mp, aes(x = dominant_metaprogram, y = pct,
                                   fill = pathologist_label)) +
  geom_col(color = "black", linewidth = 0.3) +
  scale_fill_manual(values = annot_colors, name = "Tissue Type") +
  scale_x_discrete(labels = c(MP1="MP1\nEMT", MP2="MP2\nHypoxia",
                               MP3="MP3\nSquamous", MP4="MP4\nFibroblast",
                               MP5="MP5\nKeratinocyte")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(title = "Tissue Type Composition Within Each Dominant Metaprogram",
       subtitle = "What tissue types do spots in each MP belong to?",
       x = "Dominant Metaprogram", y = "% Spots") +
  theme_pub(angle_x = 0)

ggsave(file.path(OUTPUT_DIR, "Tissue_per_DominantMP.pdf"),
       p4, width = 10, height = 6, dpi = 300)
cat("  Tissue_per_DominantMP.pdf\n")

# ==============================================================================
# FIGURE 5: Heatmap — Mean MP Score × Tissue Type × Hypoxia Status
# ==============================================================================

cat("Figure 5: Heatmap mean MP score × tissue × hypoxia...\n")

mp_heatmap_data <- mp_long %>%
  filter(!is.na(hypoxia_status)) %>%
  group_by(pathologist_label, mp, hypoxia_status) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  mutate(tissue_hypoxia = paste0(pathologist_label, "\n(", hypoxia_status, ")")) %>%
  select(tissue_hypoxia, mp, mean_score) %>%
  pivot_wider(names_from = mp, values_from = mean_score, values_fill = NA)

mat <- as.matrix(mp_heatmap_data[, -1])
rownames(mat) <- mp_heatmap_data$tissue_hypoxia

# Scale per MP (column) for better visual contrast
mat_scaled <- scale(mat)
mat_scaled[is.nan(mat_scaled)] <- 0

pdf(file.path(OUTPUT_DIR, "MP_Heatmap_Tissue_Hypoxia.pdf"), width = 8, height = max(6, nrow(mat) * 0.4))
pheatmap(
  mat_scaled,
  color       = colorRampPalette(c("#0072B2", "white", "#D55E00"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main         = "Mean MP Score by Tissue Type & Hypoxia Status\n(z-scored per MP)",
  fontsize      = 9,
  border_color  = "white",
  angle_col     = 45,
  labels_col    = c("MP1\nEMT", "MP2\nHypoxia", "MP3\nSquamous",
                    "MP4\nFibroblast", "MP5\nKeratinocyte")
)
dev.off()
cat("  MP_Heatmap_Tissue_Hypoxia.pdf\n")

# ==============================================================================
# FIGURE 6: Within-Tumor MP Shift — Hypoxic vs Normoxic
# ==============================================================================

cat("Figure 6: Tumor spots — MP shift hypoxic vs normoxic...\n")

tumor_data <- mp_long %>%
  filter(pathologist_label == "Tumor", !is.na(hypoxia_status))

if (nrow(tumor_data) < 10) {
  cat("  Insufficient tumor spots — skipping Figure 6\n")
} else {
  # Compute mean per MP × hypoxia per sample, then take difference
  tumor_shift <- tumor_data %>%
    group_by(sample, cell_line, mp, hypoxia_status) %>%
    summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = hypoxia_status, values_from = mean_score) %>%
    filter(!is.na(Hypoxic), !is.na(Normoxic)) %>%
    mutate(shift = Hypoxic - Normoxic)

  p6a <- ggplot(tumor_data, aes(x = mp, y = score, fill = hypoxia_status)) +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3,
                 position = position_dodge(width = 0.75), linewidth = 0.3) +
    scale_fill_manual(values = hypoxia_colors, name = "Region") +
    scale_x_discrete(labels = mp_labels) +
    labs(title = "Tumor Spots: MP Scores Hypoxic vs Normoxic",
         x = "Metaprogram", y = "MP Score") +
    theme_pub(angle_x = 0)

  p6b <- ggplot(tumor_shift, aes(x = mp, y = shift, fill = mp)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3, linewidth = 0.3) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
    scale_fill_manual(values = mp_colors, guide = "none") +
    scale_x_discrete(labels = mp_labels) +
    labs(title = "Hypoxia-Induced MP Shift in Tumor Spots",
         subtitle = "Per-sample: mean(Hypoxic) − mean(Normoxic)",
         x = "Metaprogram", y = "Δ Score (Hyp − Nor)") +
    theme_pub(angle_x = 0)

  combined6 <- p6a + p6b + plot_annotation(
    title = "Tumor: Metaprogram Shifts Under Hypoxia",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

  ggsave(file.path(OUTPUT_DIR, "Tumor_MP_Hypoxia_Shift.pdf"),
         combined6, width = 14, height = 6, dpi = 300)
  cat("  Tumor_MP_Hypoxia_Shift.pdf\n")
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n", strrep("═", 55), "\n", sep = "")
cat("DONE — Annotation × MP × Hypoxia plots\n\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("Files generated:\n")
cat("  1. MP_Score_by_Tissue_Type.pdf\n")
cat("  2. Hypoxia_by_Tissue_Type.pdf\n")
cat("  3. Dominant_MP_per_Tissue.pdf\n")
cat("  4. Tissue_per_DominantMP.pdf\n")
cat("  5. MP_Heatmap_Tissue_Hypoxia.pdf\n")
cat("  6. Tumor_MP_Hypoxia_Shift.pdf\n\n")
