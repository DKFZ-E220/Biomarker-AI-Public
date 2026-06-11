#!/usr/bin/env Rscript

# Resolve script directory
script_dir <- tryCatch({
  frames <- sys.frames(); ofile <- NULL
  for (f in rev(frames)) { if (!is.null(f$ofile)) { ofile <- f$ofile; break } }
  if (!is.null(ofile) && nchar(ofile) > 0) dirname(normalizePath(ofile)) else getwd()
}, error = function(e) getwd())
setwd(script_dir)

# analysis_subtypes_dir = parent of 5_BULK_SPATIAL_INTEGRATION
analysis_subtypes_dir <- dirname(script_dir)

# ==============================================================================
# BULK SUBTYPE × SPATIAL METAPROGRAM — Visualization Suite
# ==============================================================================
# 1 figure + 1 summary table:
#   Fig 1 — Stacked bar: MP composition per sample, annotated by bulk subtype
#   Summary_MP_BySubtype.csv — mean ± SD per subtype
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)

# ==============================================================================
# PATHS
# ==============================================================================

MERGED_CSV        <- file.path(script_dir, "bulk_subtype_spatial_metaprogram_merged.csv")
BASE_DIR          <- dirname(dirname(analysis_subtypes_dir))
OUTPUT_DIR        <- file.path(BASE_DIR, "Writing", "BulkSubtype_Spatial_Integration")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# COLOUR PALETTES
# ==============================================================================

# Wong palette for MPs
mp_colors <- c(
  MP1 = "#0072B2",   # Blue       — EMT/Mesenchymal
  MP2 = "#E69F00",   # Orange     — Hypoxia Response
  MP3 = "#009E73",   # Green      — Squamous Differentiation
  MP4 = "#CC79A7",   # Mauve      — ECM/Stromal
  MP5 = "#D55E00"    # Vermillion — Keratinocyte/Basal
)

mp_labels_full <- c(
  MP1 = "MP1\nEMT/Mesenchymal",
  MP2 = "MP2\nHypoxia Response",
  MP3 = "MP3\nSquamous Diff.",
  MP4 = "MP4\nECM/Stromal",
  MP5 = "MP5\nKeratinocyte/Basal"
)

# Subtype colours
subtype_colors <- c(
  Heterogenous = "#999999",
  Basal        = "#56B4E9",
  Mesenchymal  = "#D55E00",
  Atypical     = "#F0E442"
)

# Publication theme
theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title       = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = base_size,     hjust = 0.5, color = "grey40"),
      axis.title       = element_text(size = base_size,     face = "bold"),
      axis.text        = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size,     face = "bold"),
      legend.text      = element_text(size = base_size - 1),
      panel.grid.minor = element_blank(),
      strip.text       = element_text(size = base_size,     face = "bold"),
      plot.margin      = unit(c(0.4, 0.4, 0.4, 0.4), "cm")
    )
}

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("   BULK SUBTYPE × SPATIAL METAPROGRAM — VISUALIZATION SUITE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Loading data...\n")

merged <- read_csv(MERGED_CSV, show_col_types = FALSE)

# Keep only samples with spatial data
merged_sp <- merged %>% filter(!is.na(Total_Spots))

cat("Bulk samples total:", nrow(merged), "\n")
cat("Samples with spatial data:", nrow(merged_sp), "\n\n")

# ==============================================================================
# FIG 1 — Stacked bar: MP composition per sample, ordered + annotated by subtype
# ==============================================================================

cat("Fig 1: MP composition stacked bar...\n")

# Reshape to long format
mp_long <- merged_sp %>%
  select(Sample, CellLine, Subtype, MP1_pct, MP2_pct, MP3_pct, MP4_pct, MP5_pct) %>%
  pivot_longer(cols = c(MP1_pct, MP2_pct, MP3_pct, MP4_pct, MP5_pct),
               names_to = "MP", values_to = "Pct") %>%
  mutate(
    MP = str_remove(MP, "_pct"),
    # Order samples by subtype then cell line for clean grouping
    Sample_short = gsub("\\.", "\n", Sample)
  )

# Order samples: group by subtype
sample_order <- merged_sp %>%
  arrange(Subtype, CellLine, Sample) %>%
  pull(Sample)

mp_long <- mp_long %>%
  mutate(Sample = factor(Sample, levels = sample_order))

# Subtype annotation strip (top bar)
subtype_df <- merged_sp %>%
  select(Sample, Subtype) %>%
  mutate(Sample = factor(Sample, levels = sample_order),
         y = 101, yend = 105)

fig1 <- ggplot(mp_long, aes(x = Sample, y = Pct, fill = MP)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.2) +
  geom_tile(data = subtype_df,
            aes(x = Sample, y = 103, fill = NULL, colour = Subtype),
            height = 5, width = 0.85, linewidth = 1.2,
            inherit.aes = FALSE, show.legend = FALSE) +
  geom_text(data = subtype_df %>% group_by(Subtype) %>%
              summarise(Sample = Sample[ceiling(n()/2)], .groups = "drop"),
            aes(x = Sample, y = 108, label = Subtype),
            size = 3, fontface = "bold", inherit.aes = FALSE) +
  scale_fill_manual(values = mp_colors,
                    labels = c("MP1 - EMT/Mesenchymal",
                               "MP2 - Hypoxia Response",
                               "MP3 - Squamous Diff.",
                               "MP4 - ECM/Stromal",
                               "MP5 - Keratinocyte/Basal"),
                    name = "Metaprogram") +
  scale_colour_manual(values = subtype_colors) +
  scale_y_continuous(limits = c(0, 112), breaks = c(0, 25, 50, 75, 100),
                     labels = c("0", "25", "50", "75", "100%")) +
  scale_x_discrete(labels = function(x) gsub("\\.", "\n", x)) +
  labs(title = "Metaprogram Composition per Sample",
       subtitle = "Samples grouped and annotated by bulk molecular subtype",
       x = NULL, y = "% Dominant Spots") +
  theme_pub() +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        legend.position = "right")

ggsave(file.path(OUTPUT_DIR, "Fig1_MP_Composition_BySubtype.pdf"),
       fig1, width = 14, height = 6, dpi = 300)
cat("  Saved: Fig1_MP_Composition_BySubtype.pdf\n")

# ==============================================================================
# BONUS: Combined summary stats table
# ==============================================================================

cat("Writing summary stats table...\n")

summary_tbl <- merged_sp %>%
  group_by(Subtype) %>%
  summarise(
    N = n(),
    MP1_mean = round(mean(MP1_pct, na.rm=TRUE),1),
    MP1_sd   = round(sd(MP1_pct,   na.rm=TRUE),1),
    MP2_mean = round(mean(MP2_pct, na.rm=TRUE),1),
    MP2_sd   = round(sd(MP2_pct,   na.rm=TRUE),1),
    MP3_mean = round(mean(MP3_pct, na.rm=TRUE),1),
    MP3_sd   = round(sd(MP3_pct,   na.rm=TRUE),1),
    MP4_mean = round(mean(MP4_pct, na.rm=TRUE),1),
    MP4_sd   = round(sd(MP4_pct,   na.rm=TRUE),1),
    MP5_mean = round(mean(MP5_pct, na.rm=TRUE),1),
    MP5_sd   = round(sd(MP5_pct,   na.rm=TRUE),1),
    Dominant_MP = names(sort(table(DominantMP), decreasing=TRUE))[1],
    .groups = "drop"
  )

write.csv(summary_tbl,
          file.path(OUTPUT_DIR, "Summary_MP_BySubtype.csv"),
          row.names = FALSE)
cat("  Saved: Summary_MP_BySubtype.csv\n")

# ==============================================================================
# DONE
# ==============================================================================

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("ALL FIGURES COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
cat("Output:", OUTPUT_DIR, "\n\n")
cat("Files:\n")
cat("  Fig1_MP_Composition_BySubtype.pdf — stacked bar, all samples\n")
cat("  Summary_MP_BySubtype.csv          — mean ± SD table\n\n")
