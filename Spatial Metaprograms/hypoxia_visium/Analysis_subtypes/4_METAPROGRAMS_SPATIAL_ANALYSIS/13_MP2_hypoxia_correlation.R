#!/usr/bin/env Rscript

# Resolve script directory
script_dir <- tryCatch({
  frames <- sys.frames(); ofile <- NULL
  for (f in rev(frames)) { if (!is.null(f$ofile)) { ofile <- f$ofile; break } }
  if (!is.null(ofile) && nchar(ofile) > 0) dirname(normalizePath(ofile)) else getwd()
}, error = function(e) getwd())
setwd(script_dir)
analysis_subtypes_dir <- if (grepl("4_METAPROGRAMS_SPATIAL_ANALYSIS", getwd(), fixed = TRUE)) dirname(getwd()) else getwd()

# ==============================================================================
# MP2 × HYPOXIA CORRELATION ANALYSIS
# ==============================================================================
# Tests the association between continuous MP2 scores and hypoxia status at
# spot level, using Spearman correlation (per sample and overall) and
# Wilcoxon rank-sum test (hypoxic vs normoxic spots per cell line).
#
# Outputs (Writing/MP2_Hypoxia_Correlation/):
#   1. Fig_MP2_boxplot_by_cellline.pdf   — MP2 score hypoxic vs normoxic per cell line
#   2. Fig_MP2_boxplot_by_sample.pdf     — MP2 score hypoxic vs normoxic per sample
#   3. Fig_MP2_spearman_per_sample.pdf   — per-sample Spearman rho forest plot
#   4. MP2_spearman_results.csv          — per-sample rho + p-value table
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

# ==============================================================================
# CONFIG
# ==============================================================================

BASE_DIR          <- dirname(dirname(analysis_subtypes_dir))
MASTER_SPOT_TABLE <- file.path(analysis_subtypes_dir, "2_METAPROGRAMS_INTEGRATION",
                                "data", "master_spot_table.csv")
OUTPUT_DIR        <- file.path(BASE_DIR, "Writing", "MP2_Hypoxia_Correlation")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

MP2_COLOR    <- "#E69F00"   # Wong palette — MP2 Hypoxia Response
HYP_COLOR    <- "#D55E00"   # Hypoxic
NORM_COLOR   <- "#0072B2"   # Normoxic

cat(strrep("═", 63), "\n")
cat("   MP2 × HYPOXIA CORRELATION ANALYSIS\n")
cat(strrep("═", 63), "\n\n")

# ==============================================================================
# LOAD & FILTER
# ==============================================================================

master <- read_csv(MASTER_SPOT_TABLE, show_col_types = FALSE) %>%
  mutate(hypoxia_label = as.integer(hypoxia_label))

# Keep only spots with known hypoxia status and valid MP2 score
df <- master %>%
  filter(!is.na(hypoxia_label), !is.na(MP2)) %>%
  mutate(
    hypoxia_status = factor(
      ifelse(hypoxia_label == 1, "Hypoxic", "Normoxic"),
      levels = c("Normoxic", "Hypoxic")
    )
  )

# Only samples that have both hypoxic and normoxic spots
valid_samples <- df %>%
  group_by(sample) %>%
  summarise(n_hyp  = sum(hypoxia_label == 1),
            n_norm = sum(hypoxia_label == 0), .groups = "drop") %>%
  filter(n_hyp > 0, n_norm > 0) %>%
  pull(sample)

df <- df %>% filter(sample %in% valid_samples)

cat("Spots included:", nrow(df), "\n")
cat("Samples:", length(unique(df$sample)), "\n")
cat("Cell lines:", paste(sort(unique(df$cell_line)), collapse = ", "), "\n\n")

# ==============================================================================
# FIGURE 1: Boxplot — MP2 score, Hypoxic vs Normoxic, faceted by cell line
# ==============================================================================

cat("Figure 1: MP2 boxplot by cell line...\n")

# Wilcoxon p-value per cell line
wilcox_cl <- df %>%
  group_by(cell_line) %>%
  summarise(
    p_val = wilcox.test(MP2[hypoxia_status == "Hypoxic"],
                        MP2[hypoxia_status == "Normoxic"],
                        exact = FALSE)$p.value,
    y_pos = max(MP2, na.rm = TRUE) * 1.05,
    .groups = "drop"
  ) %>%
  mutate(
    sig = case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01  ~ "**",
      p_val < 0.05  ~ "*",
      p_val < 0.1   ~ "†",
      TRUE          ~ "ns"
    ),
    label = paste0(sig, "\np=", formatC(p_val, format = "e", digits = 1))
  )

p1 <- ggplot(df, aes(x = hypoxia_status, y = MP2, fill = hypoxia_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.55, linewidth = 0.4) +
  geom_text(data = wilcox_cl,
            aes(x = 1.5, y = y_pos, label = label),
            inherit.aes = FALSE, size = 2.8, fontface = "bold", vjust = 0) +
  scale_fill_manual(values = c("Normoxic" = NORM_COLOR, "Hypoxic" = HYP_COLOR),
                    guide = "none") +
  facet_wrap(~ cell_line, nrow = 1) +
  labs(
    title    = "MP2 (Hypoxia Response) Score: Hypoxic vs Normoxic Spots",
    subtitle = "Wilcoxon rank-sum test per cell line",
    x        = NULL,
    y        = "MP2 Score"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9,  hjust = 0.5, colour = "grey40"),
    strip.text    = element_text(size = 10, face = "bold"),
    axis.text.x   = element_text(size = 9),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(OUTPUT_DIR, "Fig_MP2_boxplot_by_cellline.pdf"),
       p1, width = 14, height = 6, dpi = 300)
cat("  Fig_MP2_boxplot_by_cellline.pdf\n")

# ==============================================================================
# FIGURE 2: Boxplot — MP2 score per sample (faceted)
# ==============================================================================

cat("Figure 2: MP2 boxplot by sample...\n")

wilcox_samp <- df %>%
  group_by(sample, cell_line) %>%
  summarise(
    p_val = tryCatch(
      wilcox.test(MP2[hypoxia_status == "Hypoxic"],
                  MP2[hypoxia_status == "Normoxic"],
                  exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    y_pos = max(MP2, na.rm = TRUE) * 1.05,
    .groups = "drop"
  ) %>%
  mutate(
    sig = case_when(
      is.na(p_val)  ~ "",
      p_val < 0.001 ~ "***",
      p_val < 0.01  ~ "**",
      p_val < 0.05  ~ "*",
      p_val < 0.1   ~ "†",
      TRUE          ~ "ns"
    ),
    sample_short = gsub("_SPT.*", "", sample)
  )

df <- df %>%
  mutate(sample_short = gsub("_SPT.*", "", sample))

p2 <- ggplot(df, aes(x = hypoxia_status, y = MP2, fill = hypoxia_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.55, linewidth = 0.3) +
  geom_text(data = wilcox_samp,
            aes(x = 1.5, y = y_pos, label = sig),
            inherit.aes = FALSE, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("Normoxic" = NORM_COLOR, "Hypoxic" = HYP_COLOR),
                    name = "Status") +
  facet_wrap(~ sample_short, ncol = 5) +
  labs(
    title    = "MP2 Score per Sample: Hypoxic vs Normoxic Spots",
    subtitle = "Wilcoxon rank-sum test; * p<0.05, ** p<0.01, *** p<0.001",
    x        = NULL,
    y        = "MP2 Score"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9,  hjust = 0.5, colour = "grey40"),
    strip.text    = element_text(size = 7,  face = "bold"),
    axis.text.x   = element_text(size = 7),
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(OUTPUT_DIR, "Fig_MP2_boxplot_by_sample.pdf"),
       p2, width = 16, height = ceiling(length(valid_samples) / 5) * 4 + 3,
       dpi = 300, limitsize = FALSE)
cat("  Fig_MP2_boxplot_by_sample.pdf\n")

# ==============================================================================
# SPEARMAN CORRELATION — per sample, MP2 score vs hypoxia_label (0/1)
# ==============================================================================

cat("\nFigure 3: Spearman correlation per sample...\n")

spearman_results <- df %>%
  group_by(sample, cell_line) %>%
  summarise(
    n_spots  = n(),
    rho      = cor(MP2, hypoxia_label, method = "spearman"),
    p_val    = tryCatch(
      cor.test(MP2, hypoxia_label, method = "spearman", exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups  = "drop"
  ) %>%
  mutate(
    sample_short = gsub("_SPT.*", "", sample),
    sig = case_when(
      is.na(p_val)  ~ "",
      p_val < 0.001 ~ "***",
      p_val < 0.01  ~ "**",
      p_val < 0.05  ~ "*",
      p_val < 0.1   ~ "†",
      TRUE          ~ "ns"
    ),
    fdr = p.adjust(p_val, method = "BH")
  ) %>%
  arrange(desc(rho))

cat("\n=== SPEARMAN CORRELATION: MP2 vs Hypoxia Label ===\n\n")
print(spearman_results %>%
        select(sample_short, cell_line, n_spots, rho, p_val, fdr, sig),
      n = Inf)

# Save table
write.csv(spearman_results,
          file.path(OUTPUT_DIR, "MP2_spearman_results.csv"),
          row.names = FALSE)
cat("\nMP2_spearman_results.csv saved\n")

# Overall Spearman across all spots (as a supplementary stat)
overall_cor <- cor.test(df$MP2, df$hypoxia_label,
                        method = "spearman", exact = FALSE)
cat(sprintf("\nOverall Spearman (all spots pooled): rho = %.3f, p = %.2e\n",
            overall_cor$estimate, overall_cor$p.value))
cat("(Note: use per-sample rho as primary — pooled inflates N)\n\n")

# Forest plot of per-sample rho
p3 <- ggplot(spearman_results,
             aes(x = rho, y = reorder(sample_short, rho), colour = cell_line)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_point(aes(size = n_spots), alpha = 0.85) +
  geom_text(aes(label = sig), hjust = -0.4, size = 3.5, fontface = "bold") +
  scale_colour_manual(
    values = c(FaDu  = "#0072B2", UT5   = "#E69F00",
               Cal33 = "#009E73", UT8   = "#CC79A7",
               SAS   = "#D55E00", UT45  = "#56B4E9",
               SAT   = "#F0E442"),
    name = "Cell line"
  ) +
  scale_size_continuous(name = "N spots", range = c(2, 6)) +
  labs(
    title    = "Spearman Correlation: MP2 Score vs Hypoxia (per sample)",
    subtitle = "Positive rho = higher MP2 in hypoxic spots; * p<0.05, ** p<0.01, *** p<0.001",
    x        = "Spearman ρ",
    y        = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9,  hjust = 0.5, colour = "grey40"),
    axis.text.y   = element_text(size = 8),
    legend.title  = element_text(size = 9,  face = "bold"),
    legend.text   = element_text(size = 8)
  )

ggsave(file.path(OUTPUT_DIR, "Fig_MP2_spearman_per_sample.pdf"),
       p3, width = 10, height = max(6, nrow(spearman_results) * 0.4 + 3),
       dpi = 300)
cat("Fig_MP2_spearman_per_sample.pdf\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat(strrep("═", 63), "\n")
cat("DONE\n")
cat(strrep("═", 63), "\n\n")
cat("Output directory:", OUTPUT_DIR, "\n\n")
cat("Median Spearman rho across samples:", round(median(spearman_results$rho, na.rm=TRUE), 3), "\n")
cat("Samples with rho > 0:", sum(spearman_results$rho > 0, na.rm=TRUE), "/",
    nrow(spearman_results), "\n")
cat("Samples significant (p<0.05):", sum(spearman_results$p_val < 0.05, na.rm=TRUE), "/",
    nrow(spearman_results), "\n\n")
