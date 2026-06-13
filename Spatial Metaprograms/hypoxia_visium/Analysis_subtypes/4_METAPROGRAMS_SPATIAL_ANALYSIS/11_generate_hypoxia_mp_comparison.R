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
# MP HYPOXIA vs NORMOXIA COMPARISON
# ==============================================================================
# For each sample with hypoxia data (hypoxia_label == 1 spots exist):
#   Compute % spots assigned to each MP in hypoxic vs normoxic regions
#   Shift = hypoxic% - normoxic%
#
# Outputs:
#   1. Per-sample stacked barplot (hypoxia vs normoxia composition)
#   2. Boxplot of per-sample shifts per MP
#   3. Paired t-test results table (CSV)
#   4. Combined PDF
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggrepel)

# ==============================================================================
# CONFIG
# ==============================================================================

BASE_DIR          <- dirname(dirname(analysis_subtypes_dir))
MASTER_SPOT_TABLE <- file.path(analysis_subtypes_dir, "2_METAPROGRAMS_INTEGRATION",
                                "data", "master_spot_table.csv")
OUTPUT_DIR        <- file.path(BASE_DIR, "Writing", "Hypoxia_MP_Comparison")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Wong colorblind-safe palette
mp_colors <- c(
  MP1 = "#0072B2",   # EMT/Mesenchymal
  MP2 = "#E69F00",   # Hypoxia Response
  MP3 = "#009E73",   # Squamous Differentiation
  MP4 = "#CC79A7",   # ECM/Stromal
  MP5 = "#D55E00"    # Keratinocyte/Squamous-epithelial
)
mp_labels <- c(
  MP1 = "MP1 \u2014 EMT/Mesenchymal",
  MP2 = "MP2 \u2014 Hypoxia Response",
  MP3 = "MP3 \u2014 Squamous Differentiation",
  MP4 = "MP4 \u2014 ECM/Stromal",
  MP5 = "MP5 \u2014 Keratinocyte/Squamous-epithelial"
)

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
cat("   MP HYPOXIA vs NORMOXIA COMPARISON\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")

master <- read_csv(MASTER_SPOT_TABLE, show_col_types = FALSE) %>%
  mutate(hypoxia_label = as.integer(hypoxia_label))

cat("\u2713 Total spots:", nrow(master), "\n")
cat("\u2713 Samples:", length(unique(master$sample)), "\n")

# Samples that have BOTH hypoxic and normoxic spots
samples_with_hypoxia <- master %>%
  group_by(sample) %>%
  summarise(
    n_hypoxic  = sum(hypoxia_label == 1, na.rm = TRUE),
    n_normoxic = sum(hypoxia_label == 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_hypoxic > 0, n_normoxic > 0) %>%
  pull(sample)

cat("\u2713 Samples with both hypoxic + normoxic spots:", length(samples_with_hypoxia), "\n\n")

master_hyp <- master %>%
  filter(sample %in% samples_with_hypoxia,
         !is.na(dominant_metaprogram),
         !is.na(hypoxia_label))

# ==============================================================================
# PER-SAMPLE MP COMPOSITION IN HYPOXIC vs NORMOXIC REGIONS
# ==============================================================================

cat("\ud83d\udcca Computing per-sample MP composition shifts...\n")

composition <- master_hyp %>%
  group_by(sample, cell_line, hypoxia_label) %>%
  summarise(
    total = n(),
    MP1_pct = sum(dominant_metaprogram == "MP1") / n() * 100,
    MP2_pct = sum(dominant_metaprogram == "MP2") / n() * 100,
    MP3_pct = sum(dominant_metaprogram == "MP3") / n() * 100,
    MP4_pct = sum(dominant_metaprogram == "MP4") / n() * 100,
    MP5_pct = sum(dominant_metaprogram == "MP5") / n() * 100,
    .groups = "drop"
  ) %>%
  mutate(region = ifelse(hypoxia_label == 1, "Hypoxic", "Normoxic"))

# Pivot to compute shifts
shifts_wide <- composition %>%
  select(sample, cell_line, region, MP1_pct, MP2_pct, MP3_pct, MP4_pct, MP5_pct) %>%
  pivot_wider(names_from = region,
              values_from = c(MP1_pct, MP2_pct, MP3_pct, MP4_pct, MP5_pct)) %>%
  mutate(
    MP1_shift = MP1_pct_Hypoxic - MP1_pct_Normoxic,
    MP2_shift = MP2_pct_Hypoxic - MP2_pct_Normoxic,
    MP3_shift = MP3_pct_Hypoxic - MP3_pct_Normoxic,
    MP4_shift = MP4_pct_Hypoxic - MP4_pct_Normoxic,
    MP5_shift = MP5_pct_Hypoxic - MP5_pct_Normoxic
  )

shifts_long <- shifts_wide %>%
  select(sample, cell_line, MP1_shift, MP2_shift, MP3_shift, MP4_shift, MP5_shift) %>%
  pivot_longer(cols = ends_with("_shift"),
               names_to = "MP",
               values_to = "shift_pp") %>%
  mutate(MP = str_remove(MP, "_shift"))

cat("\u2713 Shifts computed for", length(unique(shifts_wide$sample)), "samples\n\n")

# ==============================================================================
# PAIRED T-TESTS
# ==============================================================================

cat("=== PAIRED T-TEST RESULTS ===\n\n")

ttest_results <- lapply(c("MP1","MP2","MP3","MP4","MP5"), function(mp) {
  col <- paste0(mp, "_shift")
  vals <- shifts_wide[[col]]
  tt <- t.test(vals, mu = 0)
  d_val <- mean(vals, na.rm = TRUE) / sd(vals, na.rm = TRUE)  # Cohen's d vs 0
  cat(sprintf("%s: mean shift = %+.2f pp (SD %.2f), t(%d) = %.3f, p = %.4f, Cohen's d = %.3f\n",
              mp, mean(vals, na.rm=TRUE), sd(vals, na.rm=TRUE),
              tt$parameter, tt$statistic, tt$p.value, d_val))
  data.frame(
    MP           = mp,
    MP_label     = mp_labels[[mp]],
    n_samples    = sum(!is.na(vals)),
    mean_shift   = round(mean(vals, na.rm=TRUE), 2),
    sd_shift     = round(sd(vals, na.rm=TRUE), 2),
    t_stat       = round(as.numeric(tt$statistic), 3),
    df           = as.integer(tt$parameter),
    p_value      = round(tt$p.value, 4),
    cohens_d     = round(d_val, 3)
  )
}) %>% bind_rows()

cat("\n")
print(ttest_results)

# Save stats table
stats_path <- file.path(OUTPUT_DIR, "hypoxia_mp_shift_stats.csv")
write.csv(ttest_results, stats_path, row.names = FALSE)
cat("\n\u2713 Stats saved:", basename(stats_path), "\n\n")

# ==============================================================================
# CHI-SQUARE TEST: Dominant MP distribution — Hypoxic vs Normoxic
# ==============================================================================
# Tests whether the composition of dominant metaprograms differs between
# hypoxic and normoxic spots. Contingency table: rows = hypoxia status,
# cols = dominant MP (MP1–MP5).
# ==============================================================================

cat("=== CHI-SQUARE TEST: Dominant MP composition (Hypoxic vs Normoxic) ===\n\n")

# --- Overall (all spots pooled) ---
ct_overall <- table(
  ifelse(master_hyp$hypoxia_label == 1, "Hypoxic", "Normoxic"),
  master_hyp$dominant_metaprogram
)
chi_overall <- chisq.test(ct_overall)

cat(sprintf("Overall: χ²(%d) = %.2f, p = %.2e\n\n",
            chi_overall$parameter, chi_overall$statistic, chi_overall$p.value))
cat("Observed counts:\n")
print(ct_overall)
cat("\nExpected counts:\n")
print(round(chi_overall$expected, 1))
cat("\nPearson residuals (>|2| = notable deviation):\n")
print(round(chi_overall$residuals, 2))
cat("\n")

# --- Per cell line ---
cat("Per cell line:\n")
chi_cl <- lapply(sort(unique(master_hyp$cell_line)), function(cl) {
  sub <- master_hyp %>% filter(cell_line == cl)
  ct  <- table(ifelse(sub$hypoxia_label == 1, "Hypoxic", "Normoxic"),
               sub$dominant_metaprogram)
  res <- tryCatch(chisq.test(ct), error = function(e) NULL)
  if (is.null(res)) {
    data.frame(cell_line = cl, chi_stat = NA_real_,
               df = NA_integer_, p_val = NA_real_, sig = "")
  } else {
    data.frame(
      cell_line = cl,
      chi_stat  = round(as.numeric(res$statistic), 2),
      df        = as.integer(res$parameter),
      p_val     = res$p.value,
      sig       = case_when(
        res$p.value < 0.001 ~ "***",
        res$p.value < 0.01  ~ "**",
        res$p.value < 0.05  ~ "*",
        res$p.value < 0.1   ~ "†",
        TRUE                ~ "ns"
      )
    )
  }
}) %>%
  bind_rows() %>%
  mutate(fdr = p.adjust(p_val, method = "BH"))

print(chi_cl %>% select(cell_line, chi_stat, df, p_val, fdr, sig))

# Save chi-square results
chisq_out <- list(
  overall = data.frame(
    test     = "Overall (all spots pooled)",
    chi_stat = round(as.numeric(chi_overall$statistic), 2),
    df       = as.integer(chi_overall$parameter),
    p_val    = chi_overall$p.value,
    sig      = ifelse(chi_overall$p.value < 0.001, "***",
               ifelse(chi_overall$p.value < 0.01,  "**",
               ifelse(chi_overall$p.value < 0.05,  "*", "ns")))
  ),
  per_cell_line = chi_cl
)
write.csv(chisq_out$overall,
          file.path(OUTPUT_DIR, "chisq_MP_hypoxia_overall.csv"), row.names = FALSE)
write.csv(chisq_out$per_cell_line,
          file.path(OUTPUT_DIR, "chisq_MP_hypoxia_per_cellline.csv"), row.names = FALSE)
cat("\nchisq_MP_hypoxia_overall.csv\n")
cat("chisq_MP_hypoxia_per_cellline.csv\n\n")

# ==============================================================================
# FIGURE 1: Stacked bars — hypoxia vs normoxia composition per sample
# ==============================================================================

comp_plot <- composition %>%
  pivot_longer(cols = c(MP1_pct, MP2_pct, MP3_pct, MP4_pct),
               names_to = "MP", values_to = "pct") %>%
  mutate(
    MP     = str_remove(MP, "_pct"),
    region = factor(region, levels = c("Normoxic", "Hypoxic")),
    sample_short = gsub("_SPT.*", "", sample)
  )

p_bars <- ggplot(comp_plot, aes(x = region, y = pct, fill = MP)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = mp_colors, labels = mp_labels, name = "Metaprogram") +
  facet_wrap(~ sample_short, ncol = 4) +
  labs(
    title    = "Metaprogram Composition: Hypoxic vs Normoxic Regions",
    subtitle = paste0(length(samples_with_hypoxia), " samples with complete hypoxia annotations"),
    x        = NULL,
    y        = "% Spots"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9,  hjust = 0.5, colour = "grey40"),
    strip.text    = element_text(size = 7,  face = "bold"),
    axis.text.x   = element_text(size = 8),
    legend.title  = element_text(size = 9,  face = "bold"),
    legend.text   = element_text(size = 8),
    panel.grid.major.x = element_blank()
  )

# ==============================================================================
# FIGURE 2: Boxplot of shifts per MP
# ==============================================================================

# Significance labels from t-test
sig_df <- ttest_results %>%
  mutate(
    sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      p_value < 0.1   ~ "\u2020",
      TRUE            ~ "ns"
    ),
    label = sprintf("%+.1f pp\n%s", mean_shift, sig)
  )

p_box <- ggplot(shifts_long, aes(x = MP, y = shift_pp, fill = MP)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_boxplot(width = 0.5, outlier.shape = 21, outlier.size = 2, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, colour = "grey30") +
  scale_fill_manual(values = mp_colors, guide = "none") +
  scale_x_discrete(labels = mp_labels) +
  geom_text(data = sig_df,
            aes(x = MP, y = max(shifts_long$shift_pp, na.rm=TRUE) * 1.05,
                label = label),
            size = 3, fontface = "bold", vjust = 0) +
  labs(
    title    = "Metaprogram Shift: Hypoxic vs Normoxic Regions",
    subtitle = "Points = individual samples; boxes = IQR; dashed line = zero shift",
    x        = NULL,
    y        = "Shift (percentage points, Hypoxic \u2212 Normoxic)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9,  hjust = 0.5, colour = "grey40"),
    axis.text.x   = element_text(size = 9,  angle = 15, hjust = 1),
    panel.grid.major.x = element_blank()
  )

# ==============================================================================
# FIGURE 3: MP2 continuous score \u2014 Hypoxic vs Normoxic, by cell line
# ==============================================================================

cat("\ud83d\udcca Figure 3: MP2 score boxplot by cell line...\n")

HYP_COLOR  <- "#D55E00"
NORM_COLOR <- "#0072B2"

df_mp2 <- master_hyp %>%
  filter(!is.na(MP2)) %>%
  mutate(
    hypoxia_status = factor(
      ifelse(hypoxia_label == 1, "Hypoxic", "Normoxic"),
      levels = c("Normoxic", "Hypoxic")
    )
  )

wilcox_cl <- df_mp2 %>%
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
      p_val < 0.1   ~ "\u2020",
      TRUE          ~ "ns"
    ),
    label = paste0(sig, "\np=", formatC(p_val, format = "e", digits = 1))
  )

p_mp2_cl <- ggplot(df_mp2, aes(x = hypoxia_status, y = MP2, fill = hypoxia_status)) +
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
    x        = NULL, y = "MP2 Score"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle      = element_text(size = 9,  hjust = 0.5, colour = "grey40"),
    strip.text         = element_text(size = 10, face = "bold"),
    axis.text.x        = element_text(size = 9),
    panel.grid.major.x = element_blank()
  )

# ==============================================================================
# FIGURE 4: Sample-level correlation \u2014 hypoxia burden vs MP2 shift
# ==============================================================================
# Spot-level correlation is inflated (N=thousands, p always significant) and
# tests the wrong question. Instead: correlate sample-level hypoxia burden
# (% hypoxic spots per sample) with sample-level MP2 shift (hypoxic - normoxic).
# N = 17 samples \u2192 honest p-values. Tests: do samples with more hypoxia show
# a bigger MP2 shift?
# ==============================================================================

cat("\ud83d\udcca Figure 4: Per-sample MP2 score \u2014 Hypoxic vs Normoxic (paired dot plot)...\n")
# Tests whether the within-sample MP2 elevation is consistent across all samples.
# Each sample = one pair of points (mean MP2 in hypoxic vs normoxic spots), connected
# by a line. Direction of line = direction of effect. All lines going up = consistent.

sample_mp2 <- master_hyp %>%
  filter(!is.na(MP2)) %>%
  group_by(sample, cell_line, hypoxia_label) %>%
  summarise(
    mean_MP2 = mean(MP2, na.rm = TRUE),
    n_spots  = n(),
    .groups  = "drop"
  ) %>%
  mutate(
    hypoxia_status = factor(
      ifelse(hypoxia_label == 1, "Hypoxic", "Normoxic"),
      levels = c("Normoxic", "Hypoxic")
    ),
    sample_short = gsub("_SPT.*", "", sample)
  )

# Summary for annotation: how many samples have higher MP2 in hypoxic spots?
mp2_wide <- sample_mp2 %>%
  select(sample, cell_line, sample_short, hypoxia_label, mean_MP2) %>%
  pivot_wider(names_from = hypoxia_label, values_from = mean_MP2,
              names_prefix = "MP2_") %>%
  rename(MP2_hypoxic = MP2_1, MP2_normoxic = MP2_0) %>%
  mutate(direction = MP2_hypoxic > MP2_normoxic)

n_pos   <- sum(mp2_wide$direction, na.rm = TRUE)
n_total <- nrow(mp2_wide)

# Sign test: are more samples going up than chance? (binomial)
binom_p <- binom.test(n_pos, n_total, p = 0.5)$p.value
cat(sprintf("  %d/%d samples: higher MP2 in hypoxic spots (sign test p = %.4f)\n\n",
            n_pos, n_total, binom_p))

write.csv(mp2_wide,
          file.path(OUTPUT_DIR, "MP2_sample_level_consistency.csv"),
          row.names = FALSE)
cat("\u2713 MP2_sample_level_consistency.csv\n\n")

cell_line_colors <- c(FaDu = "#0072B2", UT5  = "#E69F00", Cal33 = "#009E73",
                      UT8  = "#CC79A7", SAS  = "#D55E00", UT45  = "#56B4E9",
                      SAT  = "#F0E442")

consist_label <- sprintf("%d / %d samples:\nhigher MP2 in hypoxic spots\n(sign test p = %.3f)",
                         n_pos, n_total, binom_p)

p_samplecor <- ggplot(sample_mp2,
                      aes(x = hypoxia_status, y = mean_MP2, group = sample)) +
  geom_line(aes(colour = cell_line), alpha = 0.6, linewidth = 0.7) +
  geom_point(aes(colour = cell_line, shape = cell_line), size = 2.5, alpha = 0.9) +
  annotate("text", x = 2.4, y = Inf, label = consist_label,
           hjust = 1, vjust = 1.3, size = 3.2, fontface = "bold", colour = "grey30") +
  scale_colour_manual(values = cell_line_colors, name = "Cell line") +
  scale_shape_manual(values = c(FaDu=16, UT5=17, Cal33=15, UT8=18, SAS=8, UT45=7, SAT=9),
                     name = "Cell line") +
  labs(
    title    = "MP2 Score: Hypoxic vs Normoxic Spots (per sample)",
    subtitle = "Lines connect the same sample; upward slope = higher MP2 in hypoxic spots",
    x        = NULL,
    y        = "Mean MP2 Score"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9,  hjust = 0.5, colour = "grey40"),
    legend.title  = element_text(size = 9,  face = "bold"),
    legend.text   = element_text(size = 8),
    panel.grid.major.x = element_blank()
  )

# ==============================================================================
# SAVE PDFs
# ==============================================================================

cat("\ud83d\udcc4 Saving figures...\n")

ggsave(file.path(OUTPUT_DIR, "Fig_hypoxia_mp_stacked_bars.pdf"),
       p_bars, width = 14, height = ceiling(length(samples_with_hypoxia) / 4) * 3 + 3,
       dpi = 300, limitsize = FALSE)
cat("  \u2713 Fig_hypoxia_mp_stacked_bars.pdf\n")

ggsave(file.path(OUTPUT_DIR, "Fig_hypoxia_mp_shifts_boxplot.pdf"),
       p_box, width = 10, height = 7, dpi = 300)
cat("  \u2713 Fig_hypoxia_mp_shifts_boxplot.pdf\n")

ggsave(file.path(OUTPUT_DIR, "Fig_MP2_boxplot_by_cellline.pdf"),
       p_mp2_cl, width = 14, height = 6, dpi = 300)
cat("  \u2713 Fig_MP2_boxplot_by_cellline.pdf\n")

ggsave(file.path(OUTPUT_DIR, "Fig_MP2_sample_correlation.pdf"),
       p_samplecor, width = 10, height = 7, dpi = 300)
cat("  \u2713 Fig_MP2_sample_correlation.pdf\n")

# Combined PDF
pdf(file.path(OUTPUT_DIR, "ALL_hypoxia_mp_comparison.pdf"), width = 14, height = 10)
print(p_bars)
print(p_box)
print(p_mp2_cl)
print(p_samplecor)
dev.off()
cat("  \u2713 ALL_hypoxia_mp_comparison.pdf\n\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n")
cat("\u2705 DONE\n")
cat("\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\n\n")
cat("Samples analysed:  ", length(samples_with_hypoxia), "\n")
cat("Output directory:  ", OUTPUT_DIR, "\n\n")
cat("Key findings (from t-tests):\n")
for (i in seq_len(nrow(ttest_results))) {
  r <- ttest_results[i, ]
  cat(sprintf("  %s: %+.2f pp (p = %.4f)\n", r$MP, r$mean_shift, r$p_value))
}
cat("\n")
