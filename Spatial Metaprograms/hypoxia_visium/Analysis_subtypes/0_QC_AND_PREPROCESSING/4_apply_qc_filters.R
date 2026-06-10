#!/usr/bin/env Rscript

# ==============================================================================
# APPLY PER-SAMPLE QC FILTERS
# ==============================================================================
# Filters each sample's QC Seurat object by min_UMI and min_genes thresholds
# determined from UMI distribution analysis (see UMI_distribution_stats.csv).
#
# Logic:
#   - High-quality samples (p5 > 1500): threshold 1500-2000, removes <5% of spots
#   - Decent samples (p5 > 900):        threshold 1000, removes small noise tail
#   - Samples with visible low tail:    threshold 750, cuts between two populations
#   - Problematic UT5 samples:          higher UMI + gene filter for mouse fragments
#
# Overwrites original QC RDS files in-place.
# Run BEFORE leiden.R and NMF.
# ==============================================================================

suppressPackageStartupMessages(library(Seurat))

source(here::here("config.R"))  # defines DATA_DIR, RESULTS_DIR, ANALYSIS_SUBTYPES_DIR
QC_DIR <- file.path(ANALYSIS_SUBTYPES_DIR, "0_QC_AND_PREPROCESSING/data")

# ==============================================================================
# PER-SAMPLE THRESHOLDS
# min_UMI  = minimum nCount_Spatial
# min_genes = minimum nFeature_Spatial
# ==============================================================================

thresholds <- list(
  # FaDu samples
  N150a014rep_FaDu_SPT_A1  = list(min_UMI = 1000, min_genes = 200),  # p5=986, small tail
  N150a020_FaDu_SPT_A1     = list(min_UMI = 1000, min_genes = 200),  # p5=1136, clean
  N150d320_FaDu_SPT_B1     = list(min_UMI = 1500, min_genes = 200),  # p5=1826, all high quality

  # UT5 samples — two problematic, one clean
  N153a083_UT5_C1          = list(min_UMI =  500, min_genes = 250),  # mouse fragment: gene filter critical
  N153a084_UT5_SPT_A1      = list(min_UMI = 1250, min_genes = 200),  # bimodal gap at 1000-1500
  N153b219_UT5_SPT_D1      = list(min_UMI = 1000, min_genes = 200),  # p5=1328, clean

  # Cal33 samples
  N154a037rep_Cal33_SPT_D1 = list(min_UMI = 1000, min_genes = 200),  # p5=547 but p10=1622, small tail
  N154a038_Cal33_SPT_B1    = list(min_UMI = 1500, min_genes = 200),  # p5=1833, excellent
  N154a073_Cal33_SPT_A1    = list(min_UMI = 1500, min_genes = 200),  # p5=1506, excellent

  # UT8 samples
  N155a131_UT8_SPT_B1      = list(min_UMI = 2000, min_genes = 200),  # p5=2995, exceptional
  N155b199_UT8_SPT_C1      = list(min_UMI =  750, min_genes = 200),  # low tail, gap at p25=1469
  N155b200_UT8_SPT_D1      = list(min_UMI = 1500, min_genes = 200),  # already cleaned, all pass

  # SAS samples
  N156a074_SAS_SPT_D1      = list(min_UMI = 1000, min_genes = 200),  # p5=1025, clean
  N156b140rep_SAS_SPT_B1   = list(min_UMI = 1000, min_genes = 200),  # p5=374 but small tail
  N156b181_SAS_SPT_B1      = list(min_UMI =  750, min_genes = 200),  # p10=455, gap at p25=1190

  # UT45 samples
  N157b120_UT45_SPT_C1     = list(min_UMI = 1000, min_genes = 200),  # p25=992 then jumps to p50=3904
  `N157b123_UT45-SPT_D1`   = list(min_UMI = 1000, min_genes = 200),  # p5=545, decent

  # SAT samples
  N165a002_SAT_SPT_C1      = list(min_UMI = 1500, min_genes = 200),  # p5=2123, excellent
  N165a067rep_SAT_SPT_C1   = list(min_UMI = 2000, min_genes = 200),  # p5=2388, exceptional
  N165b149rep_SAT_SPT_A1   = list(min_UMI = 1000, min_genes = 200)   # p5=653, p10=1138
)

# ==============================================================================
# APPLY FILTERS
# ==============================================================================

cat(strrep("=", 65), "\n")
cat("  APPLYING PER-SAMPLE QC FILTERS\n")
cat(strrep("=", 65), "\n\n")

results <- list()

for (sample_name in names(thresholds)) {

  rds_file <- file.path(QC_DIR, paste0(sample_name, "_qc_with_metadata.rds"))

  if (!file.exists(rds_file)) {
    cat("NOT FOUND:", sample_name, "\n\n")
    next
  }

  thr      <- thresholds[[sample_name]]
  min_umi  <- thr$min_UMI
  min_g    <- thr$min_genes

  seu <- readRDS(rds_file)
  n_before <- ncol(seu)

  seu <- seu[, seu$nCount_Spatial   >= min_umi &
               seu$nFeature_Spatial >= min_g]
  n_after  <- ncol(seu)
  n_removed <- n_before - n_after
  pct_kept  <- round(n_after / n_before * 100, 1)

  saveRDS(seu, rds_file)

  cat(sprintf("%-35s  %d → %d spots  (%.1f%% kept | removed %d)\n",
              sample_name, n_before, n_after, pct_kept, n_removed))

  results[[sample_name]] <- data.frame(
    sample    = sample_name,
    min_UMI   = min_umi,
    min_genes = min_g,
    n_before  = n_before,
    n_after   = n_after,
    pct_kept  = pct_kept,
    n_removed = n_removed
  )

  rm(seu); gc(verbose = FALSE)
}

# ==============================================================================
# SUMMARY
# ==============================================================================

summary_df <- do.call(rbind, results)
rownames(summary_df) <- NULL

out_path <- file.path(dirname(QC_DIR), "QC_inspection", "filter_summary.csv")
write.csv(summary_df, out_path, row.names = FALSE)

cat("\n", strrep("=", 65), "\n", sep = "")
cat("DONE\n")
cat(strrep("=", 65), "\n\n")
cat("Samples processed:", nrow(summary_df), "\n")
cat("Total spots removed:", sum(summary_df$n_removed), "\n")
cat("Summary saved:", out_path, "\n\n")
print(summary_df[, c("sample","n_before","n_after","pct_kept")])
