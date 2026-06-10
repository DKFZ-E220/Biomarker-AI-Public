#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

source(here::here("config.R"))  # defines DATA_DIR, RESULTS_DIR, ANALYSIS_SUBTYPES_DIR
script_dir  <- file.path(ANALYSIS_SUBTYPES_DIR, "0_QC_AND_PREPROCESSING")
qc_data_dir <- file.path(script_dir, "data")
output_dir  <- file.path(script_dir, "QC_inspection")
dir.create(output_dir, showWarnings = FALSE)

rds_files <- list.files(qc_data_dir, pattern = "_qc_with_metadata\\.rds$", full.names = TRUE)
thresholds <- c(200, 300, 500, 750, 1000, 1250, 1500, 2000, 2500, 3000)

rows <- list()

for (rds_file in rds_files) {
  sample_name <- gsub("_qc_with_metadata\\.rds$", "", basename(rds_file))
  cat("Processing:", sample_name, "\n")

  seu <- tryCatch(readRDS(rds_file), error = function(e) NULL)
  if (is.null(seu)) next

  umi   <- seu$nCount_Spatial
  genes <- seu$nFeature_Spatial

  row <- data.frame(
    sample       = sample_name,
    n_spots_total = length(umi),
    # UMI percentiles
    UMI_p5   = round(quantile(umi, 0.05)),
    UMI_p10  = round(quantile(umi, 0.10)),
    UMI_p25  = round(quantile(umi, 0.25)),
    UMI_p50  = round(quantile(umi, 0.50)),
    UMI_p75  = round(quantile(umi, 0.75)),
    UMI_p90  = round(quantile(umi, 0.90)),
    # Gene percentiles
    genes_p10 = round(quantile(genes, 0.10)),
    genes_p25 = round(quantile(genes, 0.25)),
    genes_p50 = round(quantile(genes, 0.50)),
    stringsAsFactors = FALSE
  )

  # Spots remaining at each UMI threshold
  for (thr in thresholds) {
    col <- paste0("spots_keep_UMI", thr)
    pct <- paste0("pct_keep_UMI", thr)
    row[[col]] <- sum(umi >= thr)
    row[[pct]] <- round(sum(umi >= thr) / length(umi) * 100, 1)
  }

  rows[[sample_name]] <- row
  rm(seu); gc(verbose = FALSE)
}

out <- bind_rows(rows)
out_path <- file.path(output_dir, "UMI_distribution_stats.csv")
write.csv(out, out_path, row.names = FALSE)
cat("\nSaved:", out_path, "\n")
print(out[, 1:12])
