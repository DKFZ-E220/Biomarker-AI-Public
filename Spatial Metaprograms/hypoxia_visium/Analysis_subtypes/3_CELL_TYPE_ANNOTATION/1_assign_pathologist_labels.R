#!/usr/bin/env Rscript

# ==============================================================================
# ASSIGN PATHOLOGIST LABELS TO VISIUM SPOTS
# ==============================================================================
# Area-based majority assignment: maps QuPath GeoJSON annotation polygons onto
# Visium spot capture areas from master_spot_table.csv.
# Each spot is modelled as a circle (radius = 27.5 µm / pixel_size).
# A label is assigned if ≥70% of the spot area falls within that annotation
# polygon. Spots below threshold → "Ambiguous". No overlap → "Unannotated".
# Same H&E image used in both QuPath and SpaceRanger — no registration needed.
#

#
# Output: pathologist_labels_all_samples.csv
#   columns: spot, sample, pathologist_label, majority_fraction, frac_<Label>...
# ==============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(readr)
  library(jsonlite)
  library(tidyr)
})

# ==============================================================================
# PATHS
# ==============================================================================

source(here::here("config.R"))  # defines DATA_DIR, RESULTS_DIR, ANALYSIS_SUBTYPES_DIR
BASE_DIR        <- ANALYSIS_SUBTYPES_DIR
GEOJSON_DIR     <- file.path(BASE_DIR, "3_CELL_TYPE_ANNOTATION", "geojson_exports")
SPACERANGER_DIR <- file.path(DATA_DIR, "SPT/samples_3.3.1")
MASTER_CSV      <- file.path(BASE_DIR, "2_METAPROGRAMS_INTEGRATION", "data", "master_spot_table.csv")
OUT_DIR         <- file.path(BASE_DIR, "3_CELL_TYPE_ANNOTATION", "results")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Load spot coordinates from master_spot_table (no SpaceRanger path needed)
cat("Loading spot coordinates from master_spot_table...\n")
spots_all <- read_csv(MASTER_CSV, show_col_types = FALSE) |>
  dplyr::select(barcode, sample, pxl_col_in_fullres, pxl_row_in_fullres)
cat("  Loaded", nrow(spots_all), "spots across", n_distinct(spots_all$sample), "samples\n\n")

# ==============================================================================
# SAMPLE MAPPING: N-number prefix -> full Visium sample name
# ==============================================================================

# Samples excluded from annotation (pathologist did not annotate these)
EXCLUDE_SAMPLES <- c(
  "N154a038_Cal33_SPT_B1"   # no pathologist annotation available
)

sample_map <- c(
  N150a014 = "N150a014rep_FaDu_SPT_A1",
  N150a020 = "N150a020_FaDu_SPT_A1",
  N150d320 = "N150d320_FaDu_SPT_B1",
  N153a083 = "N153a083_UT5_C1",
  N153a084 = "N153a084_UT5_SPT_A1",
  N153b219 = "N153b219_UT5_SPT_D1",
  N154a037 = "N154a037rep_Cal33_SPT_D1",
  N154a038 = "N154a038_Cal33_SPT_B1",
  N154a073 = "N154a073_Cal33_SPT_A1",
  N155a131 = "N155a131_UT8_SPT_B1",
  N155b199 = "N155b199_UT8_SPT_C1",
  N155b200 = "N155b200_UT8_SPT_D1",
  N156a074 = "N156a074_SAS_SPT_D1",
  N156b140 = "N156b140rep_SAS_SPT_B1",
  N156b181 = "N156b181_SAS_SPT_B1",
  N157b120 = "N157b120_UT45_SPT_C1",
  N157b123 = "N157b123_UT45-SPT_D1",
  N165a002 = "N165a002_SAT_SPT_C1",
  N165a067 = "N165a067rep_SAT_SPT_C1",
  N165b149 = "N165b149rep_SAT_SPT_A1"
)

# ==============================================================================
# AREA-BASED MAJORITY ASSIGNMENT PARAMETERS
# ------------------------------------------------------------------------------
# 
# Each Visium spot (55 µm diameter) is modelled as a circle. Pixel size in
# this dataset: 1 px = 0.353 µm → spot radius = 27.5 µm / 0.353 µm·px⁻¹ ≈ 77.9 px.
# For every spot we compute the intersection area with each QuPath annotation
# polygon. A label is assigned only when ≥ 70 % of the spot circle falls
# within that tissue class; otherwise the spot is labelled "Ambiguous".
# Spots with zero overlap to any annotation → "Unannotated".
# ==============================================================================

MAJORITY_THRESHOLD <- 0.70    

# Helper: read spot radius in full-res pixels from SpaceRanger scalefactors
get_spot_radius_px <- function(sample_name, spaceranger_dir) {
  scalef_path <- file.path(spaceranger_dir, sample_name,
                           "outs", "spatial", "scalefactors_json.json")
  if (!file.exists(scalef_path)) {
    stop("scalefactors_json.json not found for sample: ", sample_name,
         "\n  Expected: ", scalef_path)
  }
  scalef <- jsonlite::read_json(scalef_path)
  radius <- scalef$spot_diameter_fullres / 2
  cat("  Spot radius from SpaceRanger:", round(radius, 2), "px\n")
  radius
}

# ==============================================================================
# PROCESS EACH SAMPLE
# ==============================================================================

geojson_files <- list.files(GEOJSON_DIR, pattern = "\\.geojson$", full.names = TRUE)
all_results   <- list()

for (gjson_file in geojson_files) {
  fname  <- basename(gjson_file)
  n_code <- regmatches(fname, regexpr("N\\d+[a-z]\\d+", fname))

  if (length(n_code) == 0 || !n_code %in% names(sample_map)) {
    cat("No sample mapping for:", fname, "\n")
    next
  }

  sample_name <- sample_map[[n_code]]

  if (sample_name %in% EXCLUDE_SAMPLES) {
    cat("Skipping (no pathologist annotation):", sample_name, "\n")
    next
  }

  cat(strrep("-", 55), "\n")
  cat("Sample:", sample_name, "\n")

  # Get coordinates for this sample from master_spot_table
  pos <- spots_all |> filter(sample == sample_name)
  if (nrow(pos) == 0) {
    cat("No spots found in master_spot_table for:", sample_name, "\n\n")
    next
  }
  cat("  Spots:", nrow(pos), "\n")

  # Get per-sample spot radius from SpaceRanger scalefactors
  spot_radius_px <- get_spot_radius_px(sample_name, SPACERANGER_DIR)

  # Build sf points, then buffer to spot circles (radius = spot_radius_px)
  # QuPath uses x = col, y = row  →  same convention as pxl_col / pxl_row
  spots_pts <- st_as_sf(pos,
                        coords = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
                        crs = NA)
  spots_sf  <- st_buffer(spots_pts, dist = spot_radius_px)
  spot_area <- pi * spot_radius_px^2   # theoretical circle area in px²

  # Load annotations via jsonlite (bypasses GDAL symlink resolution issue)
  raw <- tryCatch(jsonlite::read_json(gjson_file, simplifyVector = FALSE),
                  error = function(e) NULL)
  if (is.null(raw) || length(raw) == 0) {
    cat("Failed to read GeoJSON:", fname, "\n\n")
    next
  }

  # QuPath exports as a bare JSON array of features, not a FeatureCollection
  features <- if (!is.null(raw$features)) raw$features else raw
  if (length(features) == 0) {
    cat("No features in GeoJSON:", fname, "\n\n")
    next
  }
  cat("  Annotations found:", length(features), "\n")

  # Build sf geometries directly from coordinates — no GDAL involved
  build_geom <- function(g) {
    type <- g$type
    if (type == "Polygon") {
      rings <- lapply(g$coordinates, function(ring)
        matrix(unlist(ring), ncol = 2, byrow = TRUE))
      sf::st_polygon(rings)
    } else if (type == "MultiPolygon") {
      polys <- lapply(g$coordinates, function(poly)
        lapply(poly, function(ring)
          matrix(unlist(ring), ncol = 2, byrow = TRUE)))
      sf::st_multipolygon(polys)
    } else NULL
  }

  geoms  <- list()
  labels <- character()
  for (feat in features) {
    geom <- tryCatch(build_geom(feat$geometry), error = function(e) NULL)
    if (is.null(geom)) next
    cl  <- feat$properties$classification
    lbl <- if (is.null(cl)) NA else if (is.list(cl)) cl$name else as.character(cl)
    geoms[[length(geoms) + 1]] <- geom
    labels <- c(labels, lbl)
  }

  if (length(geoms) == 0) {
    cat("No valid geometries in GeoJSON:", fname, "\n\n")
    next
  }

  annots <- sf::st_sf(label = labels, geometry = sf::st_sfc(geoms))

  # ------------------------------------------------------------------
  # AREA-BASED MAJORITY ASSIGNMENT
  # For each unique tissue label, union its polygons, intersect with
  # every spot circle, and record intersection area as a fraction of
  # the spot area.  The label whose fraction is highest and exceeds
  # MAJORITY_THRESHOLD is assigned; otherwise → "Ambiguous".
  # ------------------------------------------------------------------

  unique_labels <- sort(unique(annots$label[!is.na(annots$label)]))
  n_spots       <- nrow(spots_sf)

  # Matrix: rows = spots, cols = tissue labels  (area fractions)
  frac_mat <- matrix(0, nrow = n_spots, ncol = length(unique_labels),
                     dimnames = list(NULL, unique_labels))

  for (lbl in unique_labels) {
    polys_union <- tryCatch(
      st_union(annots[!is.na(annots$label) & annots$label == lbl, ]),
      error = function(e) NULL
    )
    if (is.null(polys_union)) next

    ix <- suppressMessages(
      st_intersects(spots_sf, polys_union, sparse = TRUE)
    )
    hit_idx <- which(lengths(ix) > 0)
    if (length(hit_idx) == 0) next

    inter_geom <- suppressMessages(
      st_intersection(spots_sf[hit_idx, ], polys_union)
    )
    inter_area <- as.numeric(st_area(inter_geom))
    frac_mat[hit_idx, lbl] <- inter_area / spot_area
  }

  # Assign label: max fraction column, only if ≥ MAJORITY_THRESHOLD
  max_frac  <- apply(frac_mat, 1, max)
  max_label <- unique_labels[apply(frac_mat, 1, which.max)]  # safe: Unannotated overrides below
  spot_labels <- dplyr::case_when(
    max_frac == 0                    ~ "Unannotated",  # no overlap at all
    max_frac >= MAJORITY_THRESHOLD   ~ max_label,      # clear winner
    TRUE                             ~ "Ambiguous"     # overlaps but no dominant label
  )

  label_counts <- table(spot_labels)
  cat("  Label distribution (area-based ≥70%):\n")
  for (l in names(label_counts)) cat(sprintf("    %-30s %d\n", l, label_counts[[l]]))

  # Build result with fraction columns for each tissue class
  frac_df <- as.data.frame(frac_mat)
  colnames(frac_df) <- paste0("frac_", gsub(" ", "_", colnames(frac_df)))

  all_results[[sample_name]] <- cbind(
    data.frame(
      spot                = pos$barcode,
      sample              = sample_name,
      pxl_col_in_fullres  = pos$pxl_col_in_fullres,
      pxl_row_in_fullres  = pos$pxl_row_in_fullres,
      pathologist_label   = spot_labels,
      majority_fraction   = max_frac,
      stringsAsFactors    = FALSE
    ),
    frac_df
  )
  cat("\n")
}

# ==============================================================================
# SAVE
# ==============================================================================

if (length(all_results) == 0) {
  stop("No samples processed. Check:\n",
       "  1. GEOJSON_DIR exists and has .geojson files\n",
       "  2. Sample names in geojson filenames match sample_map keys\n",
       "  3. master_spot_table.csv has matching sample names\n",
       "  GEOJSON_DIR: ", GEOJSON_DIR, "\n",
       "  master_spot_table: ", MASTER_CSV)
}

final_df <- dplyr::bind_rows(all_results)  # fills missing frac_ cols with NA across samples

out_path <- file.path(OUT_DIR, "pathologist_labels_all_samples.csv")
write.csv(final_df, out_path, row.names = FALSE)

cat(strrep("=", 55), "\n")
cat("Done. Total spots labelled:", nrow(final_df), "\n")
cat("Samples processed:", length(all_results), "\n")
cat("Output:", out_path, "\n\n")
print(table(final_df$pathologist_label))
