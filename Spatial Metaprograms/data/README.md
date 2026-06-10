# Data (not included in this repository)

The raw inputs are **not** committed to git. Place them under this folder (or
point `DATA_DIR` in `config.R` to wherever they live, e.g. via the
`SPATIAL_DATA_DIR` environment variable).

Expected layout under `DATA_DIR`:

```
DATA_DIR/
├── rna_sequencing/
│   └── view-by-pid/
│       └── <SAMPLE>/.../featureCounts/<SAMPLE>.fpkm_tpm.featureCounts.tsv   # bulk RNA-seq counts
└── SPT/
    ├── samples_3.3.1/<SAMPLE>/...                                           # Visium SpaceRanger outputs
    └── Hypoxia_registration_SPT_results/hypoxia_spots_updated/...           # registered hypoxia spot calls
```

All data here is preclinical (HNSCC cell-line xenograft material); no patient /
clinical data is contained in this repository.
