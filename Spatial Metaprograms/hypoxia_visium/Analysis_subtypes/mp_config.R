# ==============================================================================
# mp_config.R — Single source of truth for MP labels, colors, and sample order
# Source this at the top of every analysis script:
#   source(file.path(analysis_subtypes_dir, "mp_config.R"))
#
# UPDATE LABELS HERE after enrichment analysis confirms biological identities.
# Everything else (colors, order) stays fixed.
# ==============================================================================

# ------------------------------------------------------------------------------
# MP biological labels
# Provisional assignments based on top genes from metaprograms_final.rds:
#   MP1: VIM, AXL, ANO1, F3, EMP3        → EMT / Mesenchymal
#   MP2: DDIT4, ADM, SERPINB2, NDRG1     → Hypoxia Response
#   MP3: ECM1, CNFN, LCN2, SPRR2A        → Squamous Differentiation
#   MP4: FSTL1, FBN1, TPM1, TMSB4X       → ECM / Stromal
#   MP5: KRT4, KRTDAP, KRT13, KRT10      → Keratinocyte / Basal
# UPDATE after running 5_ENRICHMENT_ANALYSIS
# ------------------------------------------------------------------------------
mp_labels <- c(
  MP1 = "MP1 — EMT/Mesenchymal",
  MP2 = "MP2 — Hypoxia Response",
  MP3 = "MP3 — Squamous Differentiation",
  MP4 = "MP4 — ECM/Stromal",
  MP5 = "MP5 — Keratinocyte/Squamous-epithelial"
)

# Short labels for axis text / legends
mp_labels_short <- c(
  MP1 = "EMT/Mesenchymal",
  MP2 = "Hypoxia",
  MP3 = "Squamous Diff.",
  MP4 = "ECM/Stromal",
  MP5 = "Keratinocyte"
)

# ------------------------------------------------------------------------------
# Wong colorblind-safe palette — DO NOT CHANGE
# ------------------------------------------------------------------------------
wong_mp <- c(
  MP1 = "#0072B2",
  MP2 = "#E69F00",
  MP3 = "#009E73",
  MP4 = "#CC79A7",
  MP5 = "#D55E00"
)

wong_hypoxia <- c(
  Hypoxic  = "#D55E00",
  Normoxic = "#0072B2"
)

wong_cellline <- c(
  FaDu  = "#0072B2",
  SAS   = "#009E73",
  UT5   = "#D895D0",
  UT8   = "#F0E442",
  Cal33 = "#E69F00",
  SAT   = "#CC79A7"
)

# ------------------------------------------------------------------------------
# Fixed orders
# ------------------------------------------------------------------------------
mp_order        <- c("MP1", "MP2", "MP3", "MP4", "MP5")
cellline_order  <- c("FaDu", "SAS", "UT5", "UT8", "Cal33", "SAT")
hypoxia_order   <- c("Hypoxic", "Normoxic")
