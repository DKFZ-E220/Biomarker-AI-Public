# data/

Raw input data is **not** included in this repository. Clinical and molecular
data for HNSCC patients are under controlled access.

Place input files into the corresponding sub-folders listed below. Paths are
defined in `config.R` at the repository root.

## Expected sub-folder layout

```
data/
├── preclinical_transcriptomics/
│   ├── countall.csv
│   ├── count_untreated.csv
│   ├── count_treated.csv
│   └── Metadata_HNSCC.csv
├── preclinical_methylation/
│   ├── 36_Finalized_samples.csv
│   └── idat/                       # raw Illumina IDAT files
├── tcga_transcriptomics/
│   ├── Countsprimarytumormodified.csv
│   └── TCGA_metadata_GDC_RPB_combined.csv
├── tcga_methylation/
│   └── idat/                       # raw Illumina IDAT files for TCGA-HNSC
└── dktk/
    ├── dktk_adjuvant.csv
    ├── updated_centroid_values.xlsx
    ├── adjuvant_merged.csv
    └── tp53_status_dktk_adjuvant.csv
```

## Data access

| Cohort           | Source                                                         |
|------------------|----------------------------------------------------------------|
| TCGA-HNSC        | https://portal.gdc.cancer.gov (retrieved via `TCGAbiolinks`)   |
| Preclinical HNSC | Available from the authors upon reasonable request             |
| DKTK-ROG         | Controlled access — contact the authors / DKTK-ROG consortium  |

See the repository-level `README.md` for full citation and contact details.
