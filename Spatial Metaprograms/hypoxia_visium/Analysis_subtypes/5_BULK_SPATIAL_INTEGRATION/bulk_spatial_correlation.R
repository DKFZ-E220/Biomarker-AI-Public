# Bulk Molecular Subtype - Spatial Metaprogram Integration Analysis
# Integration of 20 countall bulk samples with 17 spatial metaprogram samples

# Resolve script directory — robust for source() and Rscript
script_dir <- tryCatch({
  frames <- sys.frames(); ofile <- NULL
  for (f in rev(frames)) { if (!is.null(f$ofile)) { ofile <- f$ofile; break } }
  if (!is.null(ofile) && nchar(ofile) > 0) dirname(normalizePath(ofile)) else getwd()
}, error = function(e) getwd())
setwd(script_dir)
# script_dir = 5_BULK_SPATIAL_INTEGRATION
analysis_subtypes_dir <- dirname(script_dir)
spatial_data_dir      <- dirname(dirname(analysis_subtypes_dir))  # Spatial data/

library(tidyverse)
library(readr)
library(readxl)

# ==============================================================================
# PATHS
# ==============================================================================

SUBTYPE_CSV       <- file.path(spatial_data_dir, "Bulk data", "molecular_subtypes_countall_20samples.csv")
MASTER_SPOT_TABLE <- file.path(analysis_subtypes_dir, "2_METAPROGRAMS_INTEGRATION", "data", "master_spot_table.csv")

cat("Script dir:   ", script_dir, "\n")
cat("Subtype CSV:  ", SUBTYPE_CSV, "\n")
cat("Master spots: ", MASTER_SPOT_TABLE, "\n\n")

# Load molecular subtype results (from countall centroid correlation)
subtype_results <- read_csv(SUBTYPE_CSV, show_col_types = FALSE)

# Load and aggregate master spot table (replaces old master_spot in-memory object)
master_spot <- read_csv(MASTER_SPOT_TABLE, show_col_types = FALSE) %>%
  mutate(SampleID_base = str_extract(sample, "N[0-9]+[a-z]+[0-9]+(?:rep)?"))

# Load spatial metaprogram aggregated data
# master_spot already loaded with metaprogram assignments
sample_mp_agg <- master_spot |>
  group_by(cell_line, SampleID_base) |>
  summarise(
    Total_Spots = n(),
    MP1_pct = sum(dominant_metaprogram == "MP1", na.rm = TRUE) / n() * 100,
    MP2_pct = sum(dominant_metaprogram == "MP2", na.rm = TRUE) / n() * 100,
    MP3_pct = sum(dominant_metaprogram == "MP3", na.rm = TRUE) / n() * 100,
    MP4_pct = sum(dominant_metaprogram == "MP4", na.rm = TRUE) / n() * 100,
    MP5_pct = sum(dominant_metaprogram == "MP5", na.rm = TRUE) / n() * 100,
    .groups = "drop"
  ) |>
  mutate(
    DominantMP = pmap_chr(
      list(MP1_pct, MP2_pct, MP3_pct, MP4_pct, MP5_pct),
      function(m1, m2, m3, m4, m5) {
        vals <- c(MP1=m1, MP2=m2, MP3=m3, MP4=m4, MP5=m5)
        names(which.max(vals))
      }
    )
  )

# Normalize identifiers: extract N####[a-z]### and pad numbers to 3 digits
sample_mp_agg <- sample_mp_agg |>
  mutate(
    SampleID_clean = str_remove(SampleID_base, "rep$"),  # strip trailing "rep" before parsing
    MatchID = paste0(
      str_extract(SampleID_clean, "N[0-9]+[a-z]"),
      sprintf("%03d", as.numeric(str_extract(SampleID_clean, "[0-9]+$")))
    )
  )

subtype_results <- subtype_results |>
  mutate(
    CellLine = str_extract(Sample, "^[^\\.]+"),   # e.g. "FaDu" from "FaDu.N150a.14rep"
    MatchID  = paste0(
      str_extract(Sample, "N[0-9]+[a-z]"),
      sprintf("%03d", as.numeric(sub(".*\\.(\\d+).*", "\\1", Sample)))
    )
  )

# Merge bulk subtypes with spatial metaprograms
merged_full <- subtype_results |>
  dplyr::left_join(
    sample_mp_agg |> dplyr::select(cell_line, MatchID, Total_Spots, MP1_pct, MP2_pct, MP3_pct, MP4_pct, MP5_pct, DominantMP),
    by = c("CellLine" = "cell_line", "MatchID" = "MatchID")
  ) |>
  dplyr::select(Sample, CellLine, Subtype, MaxCorrelation, Total_Spots, MP1_pct, MP2_pct, MP3_pct, MP4_pct, MP5_pct, DominantMP)

# Statistical Analysis
merged_with_spatial <- merged_full |> 
  filter(!is.na(Total_Spots))

cat("=== SUMMARY ===\n")
cat("Total bulk samples: ", nrow(merged_full), "\n")
cat("Samples with spatial data: ", nrow(merged_with_spatial), "\n")
cat("Samples without spatial data: ", nrow(merged_full) - nrow(merged_with_spatial), "\n\n")

# Contingency table
contingency <- table(merged_with_spatial$Subtype, merged_with_spatial$DominantMP)
cat("=== CONTINGENCY TABLE ===\n")
print(contingency)

# Chi-square test
chi_test <- chisq.test(contingency)
cat("\n=== CHI-SQUARE TEST ===\n")
cat("Chi-square statistic:", chi_test$statistic, "\n")
cat("p-value:", chi_test$p.value, "\n")
cat("Degrees of freedom:", chi_test$parameter, "\n")

# Descriptive statistics
stats_by_subtype <- merged_with_spatial |>
  group_by(Subtype) |>
  summarise(
    n = n(),
    MP1_mean = mean(MP1_pct, na.rm=TRUE),
    MP1_sd = sd(MP1_pct, na.rm=TRUE),
    MP2_mean = mean(MP2_pct, na.rm=TRUE),
    MP2_sd = sd(MP2_pct, na.rm=TRUE),
    MP3_mean = mean(MP3_pct, na.rm=TRUE),
    MP3_sd = sd(MP3_pct, na.rm=TRUE),
    MP4_mean = mean(MP4_pct, na.rm=TRUE),
    MP4_sd = sd(MP4_pct, na.rm=TRUE),
    MP5_mean = mean(MP5_pct, na.rm=TRUE),
    MP5_sd = sd(MP5_pct, na.rm=TRUE),
    .groups = "drop"
  )

cat("\n=== MEAN METAPROGRAM COMPOSITION BY BULK SUBTYPE ===\n\n")
print(stats_by_subtype)

# Save merged CSV
out_path <- file.path(script_dir, "bulk_subtype_spatial_metaprogram_merged.csv")
write.csv(merged_full, out_path, row.names=FALSE)
cat("\n\nResults saved to:", out_path, "\n")

# ==============================================================================
# SAVE STATISTICAL OUTPUT (chi-square + contingency table + subtype summary)
# ==============================================================================

stats_out <- file.path(script_dir, "bulk_spatial_statistical_results.txt")
sink(stats_out)

cat("================================================================================\n")
cat("  BULK SUBTYPE × DOMINANT SPATIAL METAPROGRAM — STATISTICAL RESULTS\n")
cat("  Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================================\n\n")

cat("=== SAMPLE SUMMARY ===\n")
cat("Total bulk samples:            ", nrow(merged_full), "\n")
cat("Samples with spatial data:     ", nrow(merged_with_spatial), "\n")
cat("Samples without spatial data:  ", nrow(merged_full) - nrow(merged_with_spatial), "\n\n")

cat("=== CONTINGENCY TABLE (Bulk Subtype × Dominant Spatial MP) ===\n")
print(contingency)

cat("\n=== EXPECTED COUNTS (under H0) ===\n")
print(round(chi_test$expected, 2))

cat("\n=== CHI-SQUARE TEST OF INDEPENDENCE ===\n")
cat("Chi-square statistic: ", round(chi_test$statistic, 1), "\n")
cat("Degrees of freedom:   ", chi_test$parameter, "\n")
cat("p-value:              ", formatC(chi_test$p.value, format="e", digits=3), "\n")
cat("Note: ", sum(chi_test$expected < 5), "of", length(chi_test$expected),
    "cells have expected count < 5 — chi-square approximation unreliable.\n")
cat("Result should be interpreted as descriptive/exploratory only.\n\n")

cat("=== MEAN METAPROGRAM COMPOSITION BY BULK SUBTYPE ===\n")
print(as.data.frame(stats_by_subtype))

sink()
cat("Statistical results saved to:", stats_out, "\n")

# Also save contingency table as CSV for records
ct_df <- as.data.frame.matrix(contingency)
ct_df$Subtype <- rownames(ct_df)
ct_df <- ct_df[, c("Subtype", setdiff(colnames(ct_df), "Subtype"))]
write.csv(ct_df, file.path(script_dir, "bulk_spatial_contingency_table.csv"), row.names=FALSE)
cat("Contingency table saved to:   ", file.path(script_dir, "bulk_spatial_contingency_table.csv"), "\n")
