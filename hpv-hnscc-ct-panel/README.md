# A cancer/testis antigen panel for HPV-positive head and neck squamous cell carcinoma

Multi-omics identification and clinical validation of a five-gene
cancer/testis (CT) panel (**RIBC2, TCAM1P, SMC1B, STAG3, SYCP2**) as a
prognostic biomarker for radiotherapy-treated HPV-positive head and neck
squamous cell carcinoma (HNSCC).

This repository accompanies the PhD dissertation *"[thesis title]"* by
Safayat Mahmud Khan, and contains the R / RMarkdown analysis code for the
preclinical and clinical cohorts described in the thesis.

---

## Study overview

The analysis integrates transcriptomic and methylomic data across three
cohorts, used for complementary purposes:

1. **Preclinical (cell lines)** — HPV-positive and HPV-negative HNSCC cell
   lines with matched RNA-seq and Illumina methylation data (treated and
   untreated), used for initial HPV-associated DEG / DMR discovery.
2. **TCGA-HNSC (full cohort)** — bulk RNA-seq *and* 450k/EPIC methylation
   for primary HNSCC tumours retrieved via `TCGAbiolinks`. The **entire
   TCGA cohort with HPV status** is used for the multi-omics derivation of
   the 5-gene CT panel (intersection of HPV-associated DEGs and DMRs).
3. **TCGA-HNSC (radiotherapy subset)** — the subset of TCGA-HNSC patients
   who received radiotherapy is then taken forward for unsupervised K-means
   clustering on the CT panel and survival analysis.
4. **DKTK-ROG (validation)** — an independent clinical cohort of adjuvant
   radio(chemo)therapy-treated HNSCC patients (Affymetrix HTA 2.0
   expression + long clinical follow-up) used to **validate** the
   prognostic value of the CT-panel clusters, including TP53 / p53-IHC and
   HPV-detection stratifications.

Flow in one line:
preclinical discovery → TCGA multi-omics panel derivation (whole cohort) →
survival on the RT-treated TCGA subset → DKTK-ROG validation.

---

## Repository layout

```
hpv-hnscc-ct-panel/
├── README.md
├── LICENSE                            # MIT
├── .gitignore
├── config.R                           # central path / constant config
├── data/
│   └── README.md                      # expected input files & data access
├── results/                           # outputs (gitignored)
│
├── 01_preclinical_transcriptomics/
│   ├── 01_deseq_hpv.Rmd               # DESeq2 HPV+ vs HPV- per cell line
│   ├── 02_gsea_loop.Rmd               # GSEA across MSigDB C1-C8, H
│   ├── 03_common_deg_tcga.Rmd         # intersect preclinical & TCGA DEGs
│   └── 04_gene_boxplots.Rmd           # per-gene statistical tests / plots
│
├── 02_preclinical_methylation/
│   ├── 01_rnbeads_hpv.R               # RnBeads pipeline — HPV comparison
│   ├── 02_rnbeads_treatment.R         # RnBeads pipeline — treatment arm
│   ├── 03_beta_values.Rmd             # beta -> M value conversion + filtering
│   ├── 04_common_methylation.Rmd      # intersect with TCGA methylation
│   └── 05_lola_enrichment.Rmd         # LOLA histone-mark enrichment
│
├── 03_tcga_transcriptomics/
│   ├── 01_download_tcgabiolinks.Rmd   # TCGA-HNSC STAR counts retrieval
│   └── 02_deseq_tcga.Rmd              # DESeq2 HPV+ vs HPV-
│
├── 04_tcga_methylation/
│   └── 01_rnbeads_tcga.R              # RnBeads pipeline — TCGA hg19
│
├── 05_tcga_kmeans_survival/
│   └── 01_kmeans_radiotherapy.Rmd     # K-means + survival on RT-treated TCGA
│
└── 06_dktk_validation/
    ├── 01_clustering_survival.Rmd     # K-means + OS/LRR/DMFS analysis
    ├── 02_deg_pathway.Rmd             # limma DEG + CAMERA pathway analysis
    └── 03_clinical_tp53_viz.Rmd       # TableOne, TP53 & HPV visualisations
```

Numeric prefixes give the recommended run order within each folder. The
folders themselves can be run largely independently once their respective
input data is placed under `data/`.

---

## Quick start

1. Clone the repo and open it in RStudio (the `.Rproj` file will use this
   folder as the working directory).
2. Install R ≥ 4.3 and the Bioconductor packages listed below.
3. Place your input data into `data/` following the layout in
   [`data/README.md`](data/README.md).
4. Edit `config.R` if you want to change the default `DATA_DIR` / `RESULTS_DIR`.
5. Open any `.Rmd` and knit it, or run chunk-by-chunk.

All scripts begin with:

```r
source(here::here("config.R"))
```

so no absolute paths are hard-coded — every input is read from `DATA_DIR`
and every output is written to `RESULTS_DIR`.

---

## Dependencies

R (≥ 4.3) with the following packages. The RnBeads pipelines additionally
require `RnBeads.hg19` / `RnBeads.hg38` annotation packages.

**CRAN**

```
here, tidyverse, dplyr, tidyr, readr, readxl, ggplot2, ggpubr, pheatmap,
RColorBrewer, tibble, FSA, survival, survminer, tableone, VennDiagram,
writexl, factoextra, NbClust, cluster
```

**Bioconductor**

```
DESeq2, limma, edgeR, ashr, biomaRt, TCGAbiolinks, SummarizedExperiment,
clusterProfiler, enrichplot, msigdbr, oligo, pd.hta.2.0, RnBeads,
RnBeads.hg19, RnBeads.hg38, wateRmelon, LOLA, GOstats, qvalue
```

Install with:

```r
install.packages(c("here","tidyverse","ggpubr", ...))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2","limma","TCGAbiolinks","RnBeads", ...))
```

---

## Data access

Raw input data is **not** included in this repository.

| Cohort           | How to obtain                                                |
|------------------|--------------------------------------------------------------|
| TCGA-HNSC        | Public — retrieved via `TCGAbiolinks` (see 03_tcga_*).       |
| Preclinical HNSC | Available from the authors upon reasonable request.          |
| DKTK-ROG         | Controlled access — contact the DKTK-ROG consortium / authors. |

Cluster / HPC paths that existed in the original scripts have been replaced
with relative paths driven by `config.R`. The code therefore runs on any
workstation once the data folder is populated.

---

## Reproducing the figures in the thesis

The main figures of the dissertation map onto the numbered scripts as
follows:

| Figure / Table            | Script                                             |
|---------------------------|----------------------------------------------------|
| Preclinical DEG volcano   | `01_preclinical_transcriptomics/01_deseq_hpv.Rmd`  |
| MSigDB GSEA heatmaps      | `01_preclinical_transcriptomics/02_gsea_loop.Rmd`  |
| Beta-value PCA / heatmap  | `02_preclinical_methylation/03_beta_values.Rmd`    |
| LOLA histone enrichment   | `02_preclinical_methylation/05_lola_enrichment.Rmd`|
| TCGA HPV DEG              | `03_tcga_transcriptomics/02_deseq_tcga.Rmd`        |
| TCGA K-means + KM curves  | `05_tcga_kmeans_survival/01_kmeans_radiotherapy.Rmd` |
| DKTK K-means + survival   | `06_dktk_validation/01_clustering_survival.Rmd`    |
| DKTK CAMERA pathway plot  | `06_dktk_validation/02_deg_pathway.Rmd`            |
| DKTK TP53 / HPV panels    | `06_dktk_validation/03_clinical_tp53_viz.Rmd`      |

---

## Citation


A peer-reviewed publication linking to this repository is in preparation —
the README will be updated with the DOI once available.

---

## License

Code is released under the [MIT License](LICENSE). Raw data from third-party
cohorts (TCGA, DKTK-ROG) remains governed by the data-use agreements of the
respective consortia.

---

## Contact

Questions about the code or data are best directed through the issue tracker
of this repository, or via the corresponding author of the associated
publication.
