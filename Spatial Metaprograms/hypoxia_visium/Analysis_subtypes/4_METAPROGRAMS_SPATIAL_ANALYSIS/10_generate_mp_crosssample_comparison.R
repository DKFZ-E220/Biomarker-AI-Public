#!/usr/bin/env Rscript

# Resolve script directory — robust for source() and Rscript
script_dir <- tryCatch({
  frames <- sys.frames()
  ofile <- NULL
  for (f in rev(frames)) { if (!is.null(f$ofile)) { ofile <- f$ofile; break } }
  if (!is.null(ofile) && nchar(ofile) > 0) dirname(normalizePath(ofile)) else getwd()
}, error = function(e) getwd())
setwd(script_dir)
# Normalise to Analysis_subtypes regardless of whether script_dir resolved to the R folder or its parent
analysis_subtypes_dir <- if (grepl("4_METAPROGRAMS_SPATIAL_ANALYSIS", getwd(), fixed = TRUE)) dirname(getwd()) else getwd()

# ==============================================================================
# MP CROSS-SAMPLE COMPARISON: One page per MP, all samples shown side-by-side
# ==============================================================================
# Output: Single PDF with 4 pages (one per MP)
#         Each page = all 17 samples arranged in a 5-column grid
#         Consistent global colour scale per MP for fair comparison
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

BASE_DIR          <- file.path(dirname(dirname(analysis_subtypes_dir)))
MASTER_SPOT_TABLE <- file.path(analysis_subtypes_dir, "2_METAPROGRAMS_INTEGRATION",
                                "data", "master_spot_table.csv")
OUTPUT_DIR        <- file.path(BASE_DIR, "Writing", "MP_CrossSample_Comparison")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Wong colorblind-safe palette — one colour per MP
mp_colors <- c(
  MP1 = "#0072B2",   # Blue
  MP2 = "#E69F00",   # Orange
  MP3 = "#009E73",   # Green
  MP4 = "#CC79A7",   # Mauve/Pink
  MP5 = "#D55E00"    # Vermillion
)

mp_labels <- c(
  MP1 = "MP1 — EMT/Mesenchymal",
  MP2 = "MP2 — Hypoxia Response",
  MP3 = "MP3 — Squamous Differentiation",
  MP4 = "MP4 — ECM/Stromal",
  MP5 = "MP5 — Keratinocyte/Squamous-epithelial"
)

NCOL <- 5   # columns in sample grid

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("═══════════════════════════════════════════════════════════════\n")
cat("   MP CROSS-SAMPLE COMPARISON\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Loading data...\n")
if (!file.exists(MASTER_SPOT_TABLE)) stop("Master spot table not found: ", MASTER_SPOT_TABLE)

master_data <- read_csv(MASTER_SPOT_TABLE, show_col_types = FALSE)
cat("Spots loaded:", nrow(master_data), "\n")

all_samples <- sort(unique(master_data$sample))
cat("Samples:", length(all_samples), "\n\n")

# ==============================================================================
# GLOBAL SCALES — same limits across all samples for each MP
# ==============================================================================

cat("Computing global MP scales...\n")
mp_scales <- lapply(c("MP1","MP2","MP3","MP4","MP5"), function(mp) {
  range(master_data[[mp]], na.rm = TRUE)
})
names(mp_scales) <- c("MP1","MP2","MP3","MP4","MP5")

for (mp in names(mp_scales)) {
  cat("  ", mp, ": [", round(mp_scales[[mp]][1],3), ", ", round(mp_scales[[mp]][2],3), "]\n", sep="")
}
cat("\n")

# ==============================================================================
# THEME
# ==============================================================================

theme_spatial_small <- function() {
  theme_minimal() +
    theme(
      plot.title      = element_text(size = 7, face = "bold", hjust = 0.5, margin = margin(b = 2)),
      axis.title      = element_blank(),
      axis.text       = element_blank(),
      panel.grid      = element_blank(),
      # legend kept — patchwork will collect into one shared bar
      legend.position = "right",
      legend.key.height = unit(0.35, "cm"),
      legend.key.width  = unit(0.25, "cm"),
      legend.title    = element_text(size = 7, face = "bold"),
      legend.text     = element_text(size = 6),
      plot.margin     = unit(c(2, 2, 2, 2), "pt")
    )
}

# ==============================================================================
# BUILD PAGES — one per MP
# ==============================================================================

cat("Building figures...\n")

pdf_path <- file.path(OUTPUT_DIR, "MP_CrossSample_Comparison.pdf")
pdf(pdf_path, width = 20, height = 16)

for (mp in c("MP1","MP2","MP3","MP4","MP5")) {

  cat("  Page:", mp_labels[[mp]], "...\n")

  colour_high <- mp_colors[[mp]]
  lims        <- mp_scales[[mp]]

  panel_list <- lapply(all_samples, function(s) {
    d <- master_data %>% filter(sample == s)

    # Shorten label: keep patient ID + cell line part
    short_label <- gsub("_SPT.*", "", s)   # e.g. N150a014rep_FaDu

    ggplot(d, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, colour = .data[[mp]])) +
      geom_point(size = 0.8, alpha = 0.8) +
      scale_colour_gradient(low = "white", high = colour_high,
                            limits = lims, oob = scales::squish) +
      coord_equal() +
      labs(title = short_label) +
      theme_spatial_small()
  })

  # Combine sample panels — collect identical colour guides into one shared bar
  combined <- wrap_plots(panel_list, ncol = NCOL) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title    = mp_labels[[mp]],
      subtitle = paste0("All ", length(all_samples), " samples — consistent global colour scale"),
      theme    = theme(
        plot.title    = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5, colour = "grey40"),
        plot.margin   = unit(c(0.5, 0.2, 0.2, 0.2), "cm"),
        legend.position = "right",
        legend.key.height = unit(3, "cm"),
        legend.key.width  = unit(0.5, "cm"),
        legend.title = element_text(size = 11, face = "bold"),
        legend.text  = element_text(size = 9)
      )
    )

  # Print to PDF
  print(combined)
  cat("    Done\n")
}

dev.off()

cat("\nPDF saved: MP_CrossSample_Comparison.pdf\n")
cat("Location:", OUTPUT_DIR, "\n\n")
cat("Pages:\n")
for (mp in names(mp_labels)) cat("  Page", which(names(mp_labels)==mp), "—", mp_labels[[mp]], "\n")
cat("\n")
