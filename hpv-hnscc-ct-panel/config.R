# =============================================================================
# config.R
# -----------------------------------------------------------------------------
# Central configuration file for the HPV-HNSCC CT panel analysis repository.
#
# EDIT THE PATHS BELOW ONCE, then all analysis scripts will locate their
# inputs and outputs correctly. Every .Rmd / .R file in this repository
# starts by calling:
#
#   source(here::here("config.R"))
#
# and then refers to files via file.path(DATA_DIR, "filename.csv") etc.
#
# The raw input files themselves are NOT included in this repository
# (clinical / sequencing data is under controlled access). See README.md
# for data access instructions.
# =============================================================================

# --- Root directories --------------------------------------------------------

# Folder that contains all raw input files (count matrices, metadata,
# IDATs, TP53 annotation, etc.). Place your data here before running
# any script.
DATA_DIR <- "data"

# Folder where intermediate and final result files (CSVs, PDFs, plots,
# RnBeads reports) will be written.
RESULTS_DIR <- "results"

# --- Convenience sub-folders (created on demand) -----------------------------

PRECLIN_TX_DATA   <- file.path(DATA_DIR, "preclinical_transcriptomics")
PRECLIN_METH_DATA <- file.path(DATA_DIR, "preclinical_methylation")
TCGA_TX_DATA      <- file.path(DATA_DIR, "tcga_transcriptomics")
TCGA_METH_DATA    <- file.path(DATA_DIR, "tcga_methylation")
DKTK_DATA         <- file.path(DATA_DIR, "dktk")

PRECLIN_TX_OUT   <- file.path(RESULTS_DIR, "preclinical_transcriptomics")
PRECLIN_METH_OUT <- file.path(RESULTS_DIR, "preclinical_methylation")
TCGA_TX_OUT      <- file.path(RESULTS_DIR, "tcga_transcriptomics")
TCGA_METH_OUT    <- file.path(RESULTS_DIR, "tcga_methylation")
TCGA_KM_OUT      <- file.path(RESULTS_DIR, "tcga_kmeans_survival")
DKTK_OUT         <- file.path(RESULTS_DIR, "dktk")

# --- Create output folders if missing ----------------------------------------

for (d in c(RESULTS_DIR,
            PRECLIN_TX_OUT, PRECLIN_METH_OUT,
            TCGA_TX_OUT, TCGA_METH_OUT,
            TCGA_KM_OUT, DKTK_OUT)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# --- The 5-gene CT (cancer/testis) panel -------------------------------------
# Final validated panel from the preclinical + TCGA multi-omics integration.

CT_PANEL_GENES <- c("RIBC2", "TCAM1P", "SMC1B", "STAG3", "SYCP2")
