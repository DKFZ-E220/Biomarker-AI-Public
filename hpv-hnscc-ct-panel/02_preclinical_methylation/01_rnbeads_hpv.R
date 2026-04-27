#!/usr/bin/env Rscript
# =============================================================================
# 01_rnbeads_hpv.R
# -----------------------------------------------------------------------------
# RnBeads pipeline for the preclinical HNSCC methylation cohort (EPIC v2).
# Compares HPV+ vs HPV- cell lines in two treatment arms (Untreated, Treated).
#
# Inputs   : data/preclinical_methylation/idat/        (IDAT files)
#            data/preclinical_methylation/sample_annotation.csv
# Outputs  : results/preclinical_methylation/          (RnBeads reports + .rda)
#
# HPC scheduling: this script was originally run on an LSF cluster. An
# example BSUB header is shown below — edit and uncomment as appropriate for
# your own cluster, otherwise just run as a plain Rscript.
#
# #BSUB -q long
# #BSUB -W 8:00
# #BSUB -n 6
# #BSUB -R span[hosts=1]
# #BSUB -R rusage[mem=100G]
# #BSUB -J HPV
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
})

# --- Paths -------------------------------------------------------------------

source(here::here("config.R"))

data.dir          <- PRECLIN_METH_DATA
idat.dir          <- file.path(data.dir, "idat")
sample.annotation <- file.path(data.dir, "sample_annotation.csv")

results.dir <- PRECLIN_METH_OUT
report.dir  <- results.dir

# --- RnBeads options ---------------------------------------------------------

rnb.options(
  import.idat.platform  = "probesEPICv2",
  analysis.name         = "HPV methylation",
  assembly              = "hg38",
  disk.dump.big.matrices = FALSE,
  disk.dump.bigff        = FALSE,
  logging.disk           = FALSE
)

# --- Parallel computing ------------------------------------------------------

parallel.isEnabled()
num.cores <- 6
parallel.setup(num.cores)
parallel.isEnabled()

# --- Reports directory -------------------------------------------------------

rnb.initialize.reports(report.dir)

# =============================================================================
# Data import
# =============================================================================

data.source <- c(idat.dir, sample.annotation)
result      <- rnb.run.import(
  data.source  = data.source,
  data.type    = "infinium.idat.dir",
  dir.reports  = report.dir
)
rnb.set <- result$rnb.set

summarized.regions(rnb.set)

save(rnb.set, file = file.path(report.dir, "rnb.set.DI.rda"))

# =============================================================================
# Exploratory analysis
# =============================================================================

rnb.run.exploratory(rnb.set, report.dir)
save(rnb.set, file = file.path(report.dir, "rnb.set.EA.rda"))

# =============================================================================
# Quality control
# =============================================================================

result.QC <- rnb.run.qc(rnb.set, report.dir)
save(result.QC, file = file.path(report.dir, "result.QC.rda"))
save(rnb.set,   file = file.path(report.dir, "rnb.set.QC.rda"))

# =============================================================================
# Preprocessing
# =============================================================================

rnb.options(
  assembly                          = "hg38",
  filtering.sex.chromosomes.removal = TRUE,
  identifiers.column                = "Sample",
  normalization.method              = "wm.dasen",
  normalization.background.method   = "methylumi.noob",
  normalization.plot.shifts         = TRUE,
  filtering.greedycut               = NULL,
  disk.dump.big.matrices            = FALSE,
  disk.dump.bigff                   = FALSE,
  logging.disk                      = FALSE
)

result.PC  <- rnb.run.preprocessing(rnb.set, dir.reports = report.dir)
rnb.set.PC <- result.PC$rnb.set

save(result.PC,  file = file.path(report.dir, "result.PC.rda"))
save(rnb.set.PC, file = file.path(report.dir, "rnb.set.PC.rda"))

# Set the reference level for the comparison columns.
sample.annotation.df            <- result.PC$rnb.set@pheno
sample.annotation.df$Untreated  <- factor(sample.annotation.df$Untreated,
                                          levels = c("HPV_Positive", "HPV_Negative"))
sample.annotation.df$Treated    <- factor(sample.annotation.df$Treated,
                                          levels = c("HPV_Positive", "HPV_Negative"))
result.PC$rnb.set@pheno         <- sample.annotation.df

# =============================================================================
# Differential methylation analysis
# =============================================================================

rnb.options(
  assembly                          = "hg38",
  identifiers.column                = "Sample",
  filtering.sex.chromosomes.removal = TRUE,
  differential                      = TRUE,
  differential.report.sites         = TRUE,
  differential.site.test.method     = "limma",
  differential.comparison.columns   = c("Untreated", "Treated"),
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
