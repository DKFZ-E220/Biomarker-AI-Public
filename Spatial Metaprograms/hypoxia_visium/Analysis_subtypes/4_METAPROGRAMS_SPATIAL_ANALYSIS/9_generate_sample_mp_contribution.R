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
# Figure: Sample NMF Contribution to Metaprograms
# ==============================================================================
# PURPOSE:
#   Show WHICH samples actually contributed NMF programs to the consensus MPs,
#   versus samples that only received module scores projected from other samples.
#   These are biologically different: contributors shaped the MP definition;
#   non-contributors are evaluated against a program they didn't help create.
#
# APPROACH:
#   1. Load all_programs_marker_filtered.rds — program names encode sample source
#      (format: {sample_name}.{prog_id})
#   2. Assign each program to best-matching MP via Jaccard similarity to
#      metaprograms_final.rds gene sets
#   3. Count programs contributed per sample × MP
#   4. Cross with master_spot_table.csv module scores to show:
#      — Panel A: Contribution heatmap (n_programs per sample × MP)
#      — Panel B: Module score comparison contributors vs non-contributors
#
# OUTPUTS:
#   Writing/Sample_MP_Contribution_Heatmap.pdf
#   Writing/Sample_MP_Contribution_Score_Comparison.pdf
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------------------

wong_mp <- c(
  MP1 = "#0072B2",
  MP2 = "#E69F00",
  MP3 = "#009E73",
  MP4 = "#CC79A7",
  MP5 = "#D55E00"
)

cellline_order <- c("FaDu", "SAS", "UT5", "UT8", "Cal33", "SAT")

wong_cellline <- c(
  FaDu  = "#0072B2",
  SAS   = "#009E73",
  UT5   = "#D895D0",
  UT8   = "#F0E442",
  Cal33 = "#E69F00",
  SAT   = "#CC79A7"
)

# Jaccard threshold: program must have ≥ this similarity to be assigned to an MP
JACCARD_THRESHOLD <- 0.05

BASE_DIR <- file.path(dirname(dirname(analysis_subtypes_dir)))
INTEG_DIR   <- file.path(dirname(getwd()), "data", "parallel_results", "integrated")
MASTER_SPOT <- file.path(dirname(getwd()), "2_METAPROGRAMS_INTEGRATION",
                          "data", "master_spot_table.csv")
OUTPUT_DIR  <- file.path(BASE_DIR, "Writing")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# STEP 1: Load programs and MPs
# ------------------------------------------------------------------------------

cat("Loading programs and metaprograms...\n")

programs_file <- file.path(INTEG_DIR, "all_programs_marker_filtered.rds")
if (!file.exists(programs_file)) {
  programs_file <- file.path(INTEG_DIR, "robust_programs_final.rds")
}
if (!file.exists(programs_file)) {
  programs_file <- file.path(INTEG_DIR, "robust_programs.rds")
}
cat("  Using:", basename(programs_file), "\n")

all_programs   <- readRDS(programs_file)
mp_file <- file.path(INTEG_DIR, "metaprograms_final.rds")
if (!file.exists(mp_file)) mp_file <- file.path(INTEG_DIR, "metaprogram_genes.rds")
metaprograms <- readRDS(mp_file)

cat("  Programs loaded:", length(all_programs), "\n")
cat("  Metaprograms:", names(metaprograms), "\n\n")

# ------------------------------------------------------------------------------
# STEP 2: Assign each program to its best MP via Jaccard similarity
# ------------------------------------------------------------------------------

cat("Assigning programs to MPs via Jaccard similarity...\n")

jaccard <- function(a, b) {
  i <- length(intersect(a, b))
  u <- length(union(a, b))
  if (u == 0) return(0)
  i / u
}

prog_assignments <- map_dfr(names(all_programs), function(prog_name) {
  prog_genes <- all_programs[[prog_name]]

  # Jaccard to each MP
  jac_scores <- sapply(names(metaprograms), function(mp) {
    jaccard(prog_genes, metaprograms[[mp]])
  })

  best_mp  <- names(which.max(jac_scores))
  best_jac <- max(jac_scores)

  # Extract sample name: everything before the first dot
  sample_name <- sub("\\..*", "", prog_name)

  data.frame(
    program     = prog_name,
    sample      = sample_name,
    assigned_mp = if (best_jac >= JACCARD_THRESHOLD) best_mp else NA_character_,
    jaccard     = best_jac,
    stringsAsFactors = FALSE
  )
})

cat("  Assigned:", sum(!is.na(prog_assignments$assigned_mp)),
    "of", nrow(prog_assignments), "programs (threshold:", JACCARD_THRESHOLD, ")\n\n")

# ------------------------------------------------------------------------------
# STEP 3: Count programs per sample × MP
# ------------------------------------------------------------------------------

# Pull cell line from sample name
sample_cellline <- read_csv(MASTER_SPOT, show_col_types = FALSE) %>%
  distinct(sample, cell_line)

contrib_counts <- prog_assignments %>%
  filter(!is.na(assigned_mp)) %>%
  count(sample, assigned_mp, name = "n_programs") %>%
  right_join(
    # All combinations of samples × MPs (so zeros appear)
    expand_grid(
      sample      = unique(prog_assignments$sample),
      assigned_mp = names(metaprograms)
    ),
    by = c("sample", "assigned_mp")
  ) %>%
  mutate(n_programs = replace_na(n_programs, 0)) %>%
  left_join(sample_cellline, by = "sample") %>%
  mutate(
    cell_line  = factor(cell_line, levels = cellline_order),
    is_contrib = n_programs > 0
  )

# Order samples: radioresistance, then spot count within cell line
sample_order <- read_csv(MASTER_SPOT, show_col_types = FALSE) %>%
  count(sample, cell_line, name = "n_spots") %>%
  mutate(cell_line = factor(cell_line, levels = cellline_order)) %>%
  arrange(cell_line, desc(n_spots)) %>%
  pull(sample)

contrib_counts <- contrib_counts %>%
  mutate(sample = factor(sample, levels = sample_order))

cat("=== Contribution Summary ===\n")
contrib_counts %>%
  group_by(assigned_mp) %>%
  summarise(
    n_contributors     = sum(is_contrib),
    n_non_contributors = sum(!is_contrib),
    total_programs     = sum(n_programs)
  ) %>%
  print()
cat("\n")

# ------------------------------------------------------------------------------
# PANEL A: Contribution heatmap
# ------------------------------------------------------------------------------

theme_heatmap <- function() {
  theme_minimal(base_size = 10) +
    theme(
      plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9,  hjust = 0.5, color = "grey40"),
      axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      axis.text.y   = element_text(size = 9, face = "bold"),
      axis.title    = element_text(size = 10, face = "bold"),
      panel.grid    = element_blank(),
      legend.title  = element_text(size = 9, face = "bold"),
      legend.text   = element_text(size = 8)
    )
}

p_heatmap <- ggplot(contrib_counts,
                    aes(x = sample, y = assigned_mp, fill = n_programs)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(n_programs > 0, n_programs, "")),
            size = 3, fontface = "bold", color = "white") +
  scale_fill_gradientn(
    colours = c("grey92", "#90CAF9", "#1565C0"),
    name    = "N programs\ncontributed",
    limits  = c(0, NA)
  ) +
  # Cell line annotation strip on x-axis via colour border
  # Mark contributing cells
  scale_y_discrete(
    limits = rev(names(metaprograms)),
    labels = c(
      MP1 = "MP1\nEMT/Mesen.",
      MP2 = "MP2\nHypoxia",
      MP3 = "MP3\nSquamous Diff.",
      MP4 = "MP4\nFibroblast",
      MP5 = "MP5\nKeratinocyte"
    )
  ) +
  labs(
    title    = "Sample NMF Contribution to Consensus Metaprograms",
    subtitle = paste0(
      "Numbers = NMF programs from that sample assigned to each MP\n",
      "Grey = no programs contributed (score projected from other samples)"
    ),
    x = "Sample (radioresistance order: FaDu → SAT)",
    y = "Metaprogram"
  ) +
  theme_heatmap()

ggsave(
  file.path(OUTPUT_DIR, "Sample_MP_Contribution_Heatmap.pdf"),
  p_heatmap, width = 14, height = 5, dpi = 300
)
cat("Sample_MP_Contribution_Heatmap.pdf saved\n")

# ------------------------------------------------------------------------------
# PANEL B: Module score — contributors vs non-contributors per MP
# ------------------------------------------------------------------------------

cat("Loading master_spot_table for score comparison...\n")
spots <- read_csv(MASTER_SPOT, show_col_types = FALSE)

# Contribution status per sample per MP (long form)
contrib_status <- contrib_counts %>%
  select(sample, assigned_mp, is_contrib) %>%
  rename(mp = assigned_mp)

# Melt MP scores: one row per spot × MP
scores_long <- spots %>%
  select(sample, cell_line, MP1, MP2, MP3, MP4, MP5) %>%
  pivot_longer(
    cols      = c(MP1, MP2, MP3, MP4, MP5),
    names_to  = "mp",
    values_to = "score"
  ) %>%
  left_join(contrib_status, by = c("sample", "mp")) %>%
  filter(!is.na(is_contrib)) %>%
  mutate(
    contrib_label = if_else(is_contrib, "Contributor\n(NMF programs)", "Non-contributor\n(score projected)"),
    mp_label = case_when(
      mp == "MP1" ~ "MP1 — EMT/Mesenchymal",
      mp == "MP2" ~ "MP2 — Hypoxia Response",
      mp == "MP3" ~ "MP3 — Squamous Differentiation",
      mp == "MP4" ~ "MP4 — Fibroblast/Stromal",
      mp == "MP5" ~ "MP5 — Keratinocyte/Basal"
    )
  )

# Summarise n per group for subtitle
n_summary <- scores_long %>%
  distinct(sample, mp, contrib_label) %>%
  count(mp, contrib_label)

# Per-sample means for dot overlay
sample_means <- scores_long %>%
  group_by(sample, mp, mp_label, contrib_label) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

p_scores <- ggplot(scores_long,
                   aes(x = contrib_label, y = score, fill = contrib_label)) +
  geom_violin(alpha = 0.7, trim = TRUE, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, outlier.alpha = 0.3,
               fill = "white", linewidth = 0.4) +
  # Per-sample mean dots — separate layer with its own data
  geom_jitter(
    data = sample_means,
    aes(x = contrib_label, y = mean_score),
    shape = 21, size = 2.5, fill = "grey20", color = "white",
    stroke = 0.4, width = 0.07,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Contributor\n(NMF programs)"       = "#2196F3",
      "Non-contributor\n(score projected)" = "#BDBDBD"
    ),
    guide = "none"
  ) +
  facet_wrap(~ mp_label, nrow = 1, scales = "free_y") +
  labs(
    title    = "Module Score: Contributors vs Non-Contributors",
    subtitle = paste0(
      "Dots = per-sample mean score. ",
      "Contributors shaped the MP definition; non-contributors were scored post-hoc."
    ),
    x = NULL,
    y = "Module Score (AddModuleScore)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title     = element_text(size = 12, face = "bold",  hjust = 0.5),
    plot.subtitle  = element_text(size = 9,  hjust = 0.5, color = "grey40"),
    strip.text     = element_text(size = 10, face = "bold"),
    axis.text.x    = element_text(size = 9),
    axis.title.y   = element_text(size = 10, face = "bold"),
    panel.border   = element_rect(colour = "grey80", fill = NA, linewidth = 0.4),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(OUTPUT_DIR, "Sample_MP_Contribution_Score_Comparison.pdf"),
  p_scores, width = 14, height = 6, dpi = 300
)
cat("Sample_MP_Contribution_Score_Comparison.pdf saved\n")

# ------------------------------------------------------------------------------
# BONUS: Per-MP contributor list (print to console for reference)
# ------------------------------------------------------------------------------

cat("\n=== CONTRIBUTORS PER MP ===\n")
for (mp in names(metaprograms)) {
  contributors <- contrib_counts %>%
    filter(assigned_mp == mp, is_contrib) %>%
    pull(sample) %>% as.character()
  non_contrib  <- contrib_counts %>%
    filter(assigned_mp == mp, !is_contrib) %>%
    pull(sample) %>% as.character()
  cat("\n", mp, "— Contributors (", length(contributors), "):\n")
  cat(" ", paste(contributors, collapse = ", "), "\n")
  cat(" Non-contributors (", length(non_contrib), "):\n")
  cat(" ", paste(non_contrib, collapse = ", "), "\n")
}

cat("\nAll done.\n")
cat("Outputs:\n")
cat("  Writing/Sample_MP_Contribution_Heatmap.pdf\n")
cat("  Writing/Sample_MP_Contribution_Score_Comparison.pdf\n")
