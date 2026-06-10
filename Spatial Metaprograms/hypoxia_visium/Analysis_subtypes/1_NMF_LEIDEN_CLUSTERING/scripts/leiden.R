#!/usr/bin/env Rscript
Sys.setenv(RETICULATE_PYTHON = "/usr/bin/python3")

# Load libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# Get sample file from command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript leiden.R <sample_rds_file>")
}

sample_file <- args[1]

if (!file.exists(sample_file)) {
  stop("ERROR: Sample file not found: ", sample_file)
}

sample_name <- gsub("_qc_with_metadata.rds", "", basename(sample_file))

# Reproducibility — sample-specific seed
sample_seed <- sum(utf8ToInt(sample_name)) %% .Machine$integer.max
set.seed(sample_seed)
cat("Reproducibility seed:", sample_seed, "(derived from sample name)\n")

# Directories
source(here::here("config.R"))  # defines DATA_DIR, RESULTS_DIR, ANALYSIS_SUBTYPES_DIR
output_dir <- file.path(ANALYSIS_SUBTYPES_DIR, "data/parallel_results")
leiden_dir <- file.path(output_dir, sample_name, "leiden")
dir.create(leiden_dir, recursive = TRUE, showWarnings = FALSE)

cat(strrep("=", 71), "\n")
cat("Processing sample:", sample_name, "\n")
cat("Output directory:", leiden_dir, "\n")
cat(strrep("=", 71), "\n\n")

# Load Seurat object
cat("Loading Seurat object...\n")
seu <- readRDS(sample_file)

# Normalize
cat("Normalizing data...\n")
seu <- NormalizeData(seu, assay = "Spatial", normalization.method = "LogNormalize", verbose = FALSE)

# Find variable features using ALL genes (better marker detection than 3000-gene subset)
cat("Finding variable features (all genes)...\n")
total_genes <- nrow(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = total_genes, verbose = FALSE)

# Scale all genes
cat("Scaling data...\n")
all_genes <- rownames(seu)
seu <- ScaleData(seu, features = all_genes, verbose = FALSE)

# PCA
cat("Running PCA (50 components)...\n")
seu <- RunPCA(seu, assay = "Spatial", features = all_genes, npcs = 50, verbose = FALSE)

# Use fixed 15 PCs for all samples
optimal_pcs <- 15
cat("Using", optimal_pcs, "PCs for clustering\n\n")

# Save elbow plot
cat("Saving elbow plot...\n")
pdf(file.path(leiden_dir, "elbow_plot.pdf"), width = 8, height = 6)
p_elbow <- ElbowPlot(seu, ndims = 50) +
  scale_y_continuous(trans = "log") +
  geom_vline(xintercept = optimal_pcs, color = "red", linetype = "dashed") +
  ggtitle(paste(sample_name, "- Using", optimal_pcs, "PCs")) +
  theme_bw()
print(p_elbow)
dev.off()

# Clustering
cat("Finding neighbors...\n")
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:optimal_pcs, verbose = FALSE)

cat("Running Leiden clustering...\n")
set.seed(sample_seed)
seu <- FindClusters(seu, verbose = FALSE, algorithm = 4, resolution = 0.8)

cat("Running UMAP...\n")
seu <- RunUMAP(seu, reduction = "pca", dims = 1:optimal_pcs,
               n.neighbors = 30, min.dist = 0.1, spread = 2,
               seed.use = sample_seed, verbose = FALSE)

n_clusters <- length(unique(seu$seurat_clusters))
cat("Clusters identified:", n_clusters, "\n\n")

# Save clustering plots
cat("Saving clustering plots...\n")
pdf(file.path(leiden_dir, "clustering.pdf"), width = 16, height = 6)
p1 <- DimPlot(seu, reduction = "umap", label = TRUE) +
  ggtitle(paste(sample_name, "- UMAP")) + theme_bw()
p2 <- SpatialDimPlot(seu, label = TRUE, label.size = 3) +
  ggtitle(paste(sample_name, "- Spatial"))
print(p1 + p2)
dev.off()

# FindAllMarkers — used by integrate_metaprograms.Rmd Step 1.5 for NMF program filtering
# Uses all-gene PCA clustering, giving broader marker lists than 3000-gene subset
cat("Calculating cluster marker genes for", n_clusters, "clusters...\n")
deg_df <- FindAllMarkers(
  seu,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  verbose         = FALSE
)
cat("  Marker genes found:", nrow(deg_df), "\n\n")

# Save results
# Note: seu object not saved — VisiumV2 serialization issues affect all samples
cat("Saving results...\n")

write.csv(data.frame(cell = colnames(seu), cluster = seu$seurat_clusters, sample = sample_name),
          file.path(leiden_dir, "clusters.csv"), row.names = FALSE)

write.csv(deg_df, file.path(leiden_dir, "deg_results.csv"), row.names = FALSE)

cluster_stats <- seu@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(n_cells = n(), .groups = "drop")
write.csv(cluster_stats, file.path(leiden_dir, "cluster_stats.csv"), row.names = FALSE)

metadata <- data.frame(
  sample     = sample_name,
  n_pcs_used = optimal_pcs,
  n_clusters = n_clusters,
  n_spots    = ncol(seu),
  n_genes    = nrow(seu)
)
write.csv(metadata, file.path(leiden_dir, "leiden_metadata.csv"), row.names = FALSE)

cat("\nLeiden clustering complete for", sample_name, "!\n")
cat("  Clusters:", n_clusters, "\n")
cat("  PCs used:", optimal_pcs, "\n\n")
