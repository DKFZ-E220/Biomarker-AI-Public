#!/usr/bin/env Rscript
# =============================================================================
# 01_rnbeads_tcga.R
# -----------------------------------------------------------------------------
# RnBeads pipeline (hg19) for the TCGA-HNSC methylation cohort. Compares
# HPV-positive vs HPV-negative primary tumours. Starts from a previously
# imported rnb.set.DI.rda object (produced by the initial TCGA IDAT import).
#
# Inputs   : data/tcga_methylation/idat/            (IDAT files)
#            data/tcga_methylation/rnb.set.DI.rda   (pre-imported RnBeads set)
# Outputs  : results/tcga_methylation/              (reports + .rda objects)
#
# Example LSF header (edit for your own HPC):
#
# #BSUB -q long
# #BSUB -W 8:00
# #BSUB -n 6
# #BSUB -R span[hosts=1]
# #BSUB -R rusage[mem=100G]
# #BSUB -J TCGA
# =============================================================================

suppressPackageStartupMessages({
  library(BiocParallel)
  library(GenomicFeatures)
  library(RnBeads)
  library(ggplot2)
  library(reshape2)
  library(grid)
  library(gridExtra)
  library(GOstats)
  library(LOLA)
  library(qvalue)
  library(wordcloud)
  library(hexbin)
  library(RnBeads.hg38)
  library(RnBeads.hg19)
  library(wateRmelon)
})

# --- Paths -------------------------------------------------------------------

source(here::here("config.R"))

data.dir    <- TCGA_METH_DATA
results.dir <- TCGA_METH_OUT
report.dir  <- results.dir
dir.create(report.dir, recursive = TRUE, showWarnings = FALSE)

# --- RnBeads options ---------------------------------------------------------

rnb.options(
  analysis.name          = "group methylation",
  assembly               = "hg19",
  disk.dump.big.matrices = FALSE,
  disk.dump.bigff        = FALSE,
  logging.disk           = FALSE
)

# --- Parallel computing ------------------------------------------------------

parallel.isEnabled()
num.cores <- 6
parallel.setup(num.cores)
parallel.isEnabled()

rnb.initialize.reports(report.dir)

# Load the previously saved dataset
load(file = file.path(data.dir, "rnb.set.DI.rda"))

# =============================================================================
# Preprocessing
# =============================================================================

rnb.options(
  assembly                          = "hg19",
  filtering.sex.chromosomes.removal = TRUE,
  identifiers.column                = "Sample",
  filtering.greedycut               = NULL,
  disk.dump.big.matrices            = FALSE,
  disk.dump.bigff                   = FALSE,
  logging.disk                      = FALSE
)

result.PC  <- rnb.run.preprocessing(rnb.set, dir.reports = report.dir)
rnb.set.PC <- result.PC$rnb.set

save(result.PC,  file = file.path(report.dir, "result.PC.rda"))
save(rnb.set.PC, file = file.path(report.dir, "rnb.set.PC.rda"))

# Set the reference level for the HPV Subtype factor.
sample.annotation.df          <- result.PC$rnb.set@pheno
sample.annotation.df$Subtype  <- factor(sample.annotation.df$Subtype,
                                        levels = c("HPV_Positive", "HPV_Negative"))
result.PC$rnb.set@pheno       <- sample.annotation.df

# =============================================================================
# Differential methylation analysis
# =============================================================================

rnb.options(
  assembly                          = "hg19",
  identifiers.column                = "Sample",
  filtering.sex.chromosomes.removal = TRUE,
  differential                      = TRUE,
  differential.report.sites         = TRUE,
  differential.site.test.method     = "limma",
  differential.comparison.columns   = ("Subtype"),
  disk.dump.big.matrices            = FALSE,
  disk.dump.bigff                   = FALSE,
  logging.disk                      = FALSE,
  differential.enrichment.go        = TRUE,
  differential.enrichment.lola      = TRUE
)

rnb.set.DMA <- rnb.run.differential(result.PC$rnb.set, report.dir)
save(rnb.set.DMA, file = file.path(report.dir, "rnb.set.DMA.rda"))

parallel.teardown()
parallel.isEnabled()
