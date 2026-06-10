#!/usr/bin/env Rscript
Sys.setenv(RETICULATE_PYTHON = "/usr/bin/python3")

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(NMF)
  library(dplyr)
  library(Matrix)
  library(reticulate)
})

# Get sample file from command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript hypoxia_parallel_analysis_brunet.R <sample_rds_file>")
}

sample_file <- args[1]

# ...existing code for debug, sample_name, seed, directories...
cat("Raw argument:", sample_file, "\n")
cat("File exists:", file.exists(sample_file), "\n")

if (!file.exists(sample_file)) {
  stop("ERROR: Sample file not found: ", sample_file)
}

sample_name <- gsub("_qc_with_metadata.rds", "", basename(sample_file))
sample_seed <- sum(utf8ToInt(sample_name)) %% .Machine$integer.max
set.seed(sample_seed)
cat("Sample-specific seed:", sample_seed, "\n")

source(here::here("config.R"))  # defines DATA_DIR, RESULTS_DIR, ANALYSIS_SUBTYPES_DIR
out_dir <- file.path(ANALYSIS_SUBTYPES_DIR, "data/parallel_results")
dir.create(file.path(out_dir, sample_name), recursive = TRUE, showWarnings = FALSE)

cat("Processing sample:", sample_name, "\n")
cat("Input file:", sample_file, "\n")
cat("Output directory:", file.path(out_dir, sample_name), "\n\n")

seu <- readRDS(sample_file)

# ===== NMF ANALYSIS =====
# Note: Leiden clustering is run separately via leiden.R / submit_leiden_all20.sh
# using all-gene PCA for better marker gene detection
cat("Running NMF analysis...\n")

if (DefaultAssay(seu) == "Spatial") {
  layer_count <- length(Layers(seu, assay = "Spatial"))
  if (layer_count > 1) {
    cat("Joining", layer_count, "layers...\n")
    seu[["Spatial"]] <- JoinLayers(seu[["Spatial"]])
  }
}

mat <- LayerData(seu, assay = "Spatial", layer = "data")

if (nrow(mat) == 0) {
  cat("Normalizing from counts...\n")
  seu <- NormalizeData(seu, assay = "Spatial", normalization.method = "LogNormalize")
  mat <- LayerData(seu, assay = "Spatial", layer = "data")
}

mat <- as.matrix(mat)
cat("Raw matrix: Genes =", nrow(mat), "| Spots =", ncol(mat), "\n")

mat <- mat[!grepl("^MT-|^RPL|^RPS", rownames(mat)), ]
cat("After MT/ribo filter: Genes =", nrow(mat), "\n")

gene_var <- apply(mat, 1, var)
mat <- mat[gene_var > 0, ]
cat("After variance filter: Genes =", nrow(mat), "\n")

expressed_genes <- rownames(mat)[rowMeans(mat) > 0.4]
mat <- mat[expressed_genes, ]
cat("After expression filter (mean>0.4): Genes =", nrow(mat), "\n")

mat_centered <- mat - rowMeans(mat)
mat_centered[mat_centered < 0] <- 0
mat_centered[is.na(mat_centered)] <- 0
mat_centered[is.infinite(mat_centered)] <- 0

cat("Final matrix: Genes =", nrow(mat_centered), "| Spots =", ncol(mat_centered), "\n")
cat("Matrix stats: min =", min(mat_centered), "max =", max(mat_centered), "mean =", mean(mat_centered), "\n")

genes_with_signal <- rowSums(mat_centered) > 0
if (sum(genes_with_signal) < 50) {
  cat("ERROR: Too few genes with signal\n")
  quit(status = 0)
}

mat_centered <- mat_centered[genes_with_signal, ]
cat("After removing zero genes: Genes =", nrow(mat_centered), "\n")

# Additional matrix validation - check for problematic genes/spots
cat("Checking for problematic genes and spots...\n")

# Remove genes with extremely low variance (can cause matrix issues)
gene_sd <- apply(mat_centered, 1, sd)
valid_genes <- gene_sd > 1e-10
cat("Genes with valid SD:", sum(valid_genes), "/", length(valid_genes), "\n")
mat_centered <- mat_centered[valid_genes, ]

# Remove spots with all/mostly zeros
spot_nonzero <- colSums(mat_centered > 0)
valid_spots <- spot_nonzero > (nrow(mat_centered) * 0.01)  # At least 1% of genes expressed
cat("Valid spots:", sum(valid_spots), "/", length(valid_spots), "\n")
mat_centered <- mat_centered[, valid_spots]

# Ensure matrix is strictly non-negative and numeric
# Save names before matrix() conversion — it drops them
rn <- rownames(mat_centered)
cn <- colnames(mat_centered)
mat_centered <- pmax(mat_centered, 0)
mat_centered <- matrix(as.numeric(mat_centered), nrow=nrow(mat_centered), ncol=ncol(mat_centered))
rownames(mat_centered) <- rn
colnames(mat_centered) <- cn

cat("Final cleaned matrix: Genes =", nrow(mat_centered), "| Spots =", ncol(mat_centered), "\n")
cat("Matrix class:", class(mat_centered), "| Mode:", mode(mat_centered), "\n")
cat("Any NA:", any(is.na(mat_centered)), "| Any Inf:", any(is.infinite(mat_centered)), "\n")

nmf_dir <- file.path(out_dir, sample_name, "nmf")
dir.create(nmf_dir, showWarnings = FALSE)

cat("Running NMF rank-by-rank with nrun=5, method=brunet...\n")

nmf_results <- list()
all_programs <- list()

for (k in 2:11) {
  cat("\n=== Processing rank k=", k, " ===\n")
  rank_file <- file.path(nmf_dir, paste0("nmf_rank_", k, ".rds"))
  
  if (file.exists(rank_file)) {
    cat("Rank", k, "already exists, loading...\n")
    fit_k <- readRDS(rank_file)
  } else {
    start_time <- Sys.time()
    fit_k <- tryCatch({
      nmf(x = mat_centered, rank = k, nrun = 5, method = "brunet", seed = sample_seed)
    }, error = function(e) {
      cat("ERROR:", conditionMessage(e), "\n")
      return(NULL)
    })
    
    elapsed <- difftime(Sys.time(), start_time, units="mins")
    cat("Rank", k, "completed in", round(elapsed, 1), "minutes\n")
    
    if (is.null(fit_k)) {
      next
    }
    
    saveRDS(fit_k, rank_file)
    cat("Saved:", rank_file, "\n")
  }
  
  # Extract programs - FIX for NMFfitX1 objects
  cat("Extracting programs from rank", k, "...\n")
  
  basis_mat <- tryCatch({
    # For NMFfitX1 objects, extract best fit first
    if (inherits(fit_k, "NMFfitX1") || inherits(fit_k, "NMFfitX")) {
      cat("Extracting best fit from", class(fit_k), "...\n")
      best_fit <- fit(fit_k)
      basis(best_fit)
    } else {
      basis(fit_k)
    }
  }, error = function(e) {
    cat("ERROR extracting basis:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (is.null(basis_mat) || !is.matrix(basis_mat) || ncol(basis_mat) == 0) {
    cat("WARNING: Invalid basis matrix for rank", k, "\n")
    next
  }

  # Cached rank files were saved with NULL rownames — restore from current mat_centered
  if (is.null(rownames(basis_mat)) && nrow(basis_mat) == nrow(mat_centered)) {
    rownames(basis_mat) <- rownames(mat_centered)
    cat("  Restored", nrow(basis_mat), "gene names from mat_centered\n")
  }

  cat("Extracting", ncol(basis_mat), "programs from rank", k, "...\n")

  for (i in seq_len(ncol(basis_mat))) {
    prog_name <- paste0(sample_name, ".", k, ".P", i)
    gene_weights <- basis_mat[, i]
    top_genes <- names(sort(gene_weights, decreasing = TRUE))[1:min(50, length(gene_weights))]
    all_programs[[prog_name]] <- top_genes
  }
  
  nmf_results[[as.character(k)]] <- fit_k
}

saveRDS(nmf_results, file.path(nmf_dir, "nmf_results.rds"))
saveRDS(all_programs, file.path(nmf_dir, "all_programs.rds"))

cat("\nNMF complete. Total programs:", length(all_programs), "\n")

if (length(all_programs) == 0) {
  quit(status = 1)
}

nmf_stats <- data.frame(
  sample = sample_name,
  n_programs = length(all_programs),
  n_genes = nrow(mat_centered),
  n_cells = ncol(seu)
)
write.csv(nmf_stats, file.path(nmf_dir, "nmf_stats.csv"), row.names = FALSE)

cat("Sample processing complete:", sample_name, "\n\n")
