# =============================================================================
# config.R
# -----------------------------------------------------------------------------
# Central path configuration for the Spatial Metaprograms analysis repository.
#
# EDIT DATA_DIR BELOW ONCE (or set the SPATIAL_DATA_DIR environment variable),
# then every analysis script locates its inputs and outputs correctly. Each
# .R / .Rmd file in this repository begins with:
#
#   source(here::here("config.R"))
#
# and then refers to files via file.path(DATA_DIR, ...) etc.
#
# Raw sequencing / Visium data is NOT included in this repository. See
# data/README.md for the expected inputs and their layout.
# =============================================================================

# --- Repository root ---------------------------------------------------------
# Resolved automatically (looks for the .here / .git marker at the repo root).
if (requireNamespace("here", quietly = TRUE)) {
  SPATIAL_ROOT <- here::here()
} else {
  SPATIAL_ROOT <- getwd()
}

# --- External raw data root --------------------------------------------------
# Point this at the folder that holds the raw inputs on YOUR machine:
#   - rna_sequencing/view-by-pid/...      (bulk featureCounts .tsv files)
#   - SPT/samples_3.3.1/...               (Visium SpaceRanger outputs)
#   - SPT/Hypoxia_registration_SPT_results/...
# Defaults to <repo>/data if the SPATIAL_DATA_DIR env var is not set.
DATA_DIR <- Sys.getenv("SPATIAL_DATA_DIR", unset = file.path(SPATIAL_ROOT, "data"))

# --- Repository-internal locations -------------------------------------------
# Code, intermediate tables, and figures live under here.
ANALYSIS_SUBTYPES_DIR <- file.path(SPATIAL_ROOT, "hypoxia_visium", "Analysis_subtypes")

# Where generated results / figures are written (git-ignored).
RESULTS_DIR <- file.path(SPATIAL_ROOT, "results")

# --- Create output folders if missing ----------------------------------------
for (d in c(RESULTS_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
