#!/usr/bin/env Rscript

# ==============================================================================
# VISUALIZE PATHOLOGIST ANNOTATION vs HYPOXIA — ALL SAMPLES
# ==============================================================================
# Side-by-side spatial plots per sample:
#   Left  — spots coloured by pathologist label (QuPath annotation)
#   Right — spots coloured by hypoxia label (PIMO co-registration)
# Produces one PDF per sample + one summary PDF with all samples.
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

source(here::here("config.R"))  # defines DATA_DIR, RESULTS_DIR, ANALYSIS_SUBTYPES_DIR
ANNOT_FILE  <- file.path(ANALYSIS_SUBTYPES_DIR, "3_CELL_TYPE_ANNOTATION/results/pathologist_labels_all_samples.csv")
HYPOXIA_FILE <- file.path(ANALYSIS_SUBTYPES_DIR, "2_METAPROGRAMS_INTEGRATION/data/all_samples_aligned_hypoxia.csv")
OUT_DIR     <- file.path(ANALYSIS_SUBTYPES_DIR, "3_CELL_TYPE_ANNOTATION/results/annotation_check_plots")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# COLOURS
# ==============================================================================

label_colours <- c(
  "Tumor"                   = "#E69F00",
  "Necrosis"                = "#000000",
  "Stroma"                  = "#56B4E9",
  "Keratinization"          = "#009E73",
  "Blood vessels"           = "#D55E00",
  "Inflammatory infiltrate" = "#CC79A7",
  "Immune cells"            = "#CC79A7",
  "Edema"                   = "#0072B2",
  "Muscle"                  = "#8B4513",
  "Epidermis"               = "#F0E442",
  "Bacteria"                = "#999999",
  "Ambiguous"               = "#AAAAAA",
  "Unannotated"             = "#DDDDDD"
)

hypoxia_colours <- c("1" = "#D55E00", "0" = "#0072B2", "NA" = "#DDDDDD")

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading annotation data...\n")
annot <- read.csv(ANNOT_FILE, stringsAsFactors = FALSE)

cat("Loading hypoxia data...\n")
hypoxia <- read.csv(HYPOXIA_FILE, stringsAsFactors = FALSE) %>%
  select(barcode, sample, hypoxia_label) %>%
  rename(spot = barcode)

# Merge
df <- annot %>%
  left_join(hypoxia, by = c("spot", "sample")) %>%
  mutate(
    hypoxia_label = case_when(
      hypoxia_label == 1  ~ "1",
      hypoxia_label == 0  ~ "0",
      TRUE                ~ "NA"
    ),
    pathologist_label = factor(pathologist_label, levels = names(label_colours),
                               exclude = NULL)
  )

samples <- unique(df$sample)
cat("Samples:", length(samples), "\n\n")

# ==============================================================================
# PLOT FUNCTION
# ==============================================================================

plot_sample <- function(sdf, sample_name) {
  # flip y axis to match image orientation (pixel row increases downward)
  p1 <- ggplot(sdf, aes(x = pxl_col_in_fullres, y = -pxl_row_in_fullres,
                         colour = pathologist_label)) +
    geom_point(size = 0.6) +
    scale_colour_manual(values = label_colours, na.value = "#DDDDDD",
                        name = "Pathologist") +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "right",
          legend.text = element_text(size = 7),
          plot.title = element_text(size = 8, face = "bold")) +
    ggtitle("Pathologist annotation")

  p2 <- ggplot(sdf, aes(x = pxl_col_in_fullres, y = -pxl_row_in_fullres,
                         colour = hypoxia_label)) +
    geom_point(size = 0.6) +
    scale_colour_manual(values = hypoxia_colours,
                        labels = c("1" = "Hypoxic", "0" = "Normoxic", "NA" = "No data"),
                        name = "Hypoxia") +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "right",
          legend.text = element_text(size = 7),
          plot.title = element_text(size = 8, face = "bold")) +
    ggtitle("Hypoxia (PIMO)")

  p1 + p2 +
    plot_annotation(title = sample_name,
                    theme = theme(plot.title = element_text(size = 10, face = "bold")))
}

# ==============================================================================
# PER-SAMPLE PDFs
# ==============================================================================

cat("Generating per-sample plots...\n")

for (s in samples) {
  sdf <- df %>% filter(sample == s)
  p   <- plot_sample(sdf, s)

  out_pdf <- file.path(OUT_DIR, paste0(s, "_annotation_check.pdf"))
  ggsave(out_pdf, p, width = 14, height = 6)
  cat("  Saved:", basename(out_pdf), "\n")
}

# ==============================================================================
# SUMMARY PDF — all samples on one file
# ==============================================================================

cat("\nGenerating summary PDF...\n")

pdf(file.path(OUT_DIR, "ALL_SAMPLES_annotation_check.pdf"), width = 14, height = 6)
for (s in samples) {
  sdf <- df %>% filter(sample == s)
  p   <- plot_sample(sdf, s)
  print(p)
}
dev.off()

cat("\nDone. Plots saved to:", OUT_DIR, "\n")
