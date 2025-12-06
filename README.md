# ICB Raw Data Reproduction Project

## Overview
This project reproduces the single-cell RNA-seq analysis from the paper:
**"Integrated cancer cell-specific single-cell RNA-seq datasets"** (Gondal et al.)

We downloaded raw data from GEO and processed it following the paper methodology.

## Datasets

| Dataset | GEO Accession | Cancer Type | Our Cells | Paper Cells |
|---------|---------------|-------------|-----------|-------------|
| Melanoma | GSE115978 | Skin melanoma | 2,018 | 1,945 |
| BCC | GSE123813 | Basal cell carcinoma | 3,452 | 3,500 |
| Liver | GSE125449 | HCC/iCCA | 1,577 | 1,992 |

## Analysis Pipeline

1. **Data Download** - Raw counts from GEO
2. **Seurat Object Creation** - CreateSeuratObject()
3. **Cell Type Filtering** - Keep malignant/tumor cells only
4. **DoubletFinder** - Remove doublets (10X data only)
5. **Standard Processing** - Normalize, Scale, PCA, UMAP, Clustering
6. **Differential Expression** - Pre vs Post treatment comparison
7. **Visualization** - UMAP and Volcano plots

## Directory Structure
```
ICB_raw_reproduction/
├── data/
│   ├── raw/
│   │   ├── melanoma/
│   │   ├── bcc/
│   │   └── liver/
│   └── processed/
│       ├── melanoma_final_processed.RDS
│       ├── bcc_final_processed.RDS
│       └── liver_final_processed.RDS
├── results/
│   ├── figures/
│   │   ├── melanoma_umap.png
│   │   ├── bcc_umap.png
│   │   ├── liver_umap.png
│   │   ├── volcano_melanoma.png
│   │   ├── volcano_bcc.png
│   │   └── volcano_liver.png
│   └── tables/
│       ├── melanoma_DE_post_vs_pre.csv
│       ├── bcc_DE_post_vs_pre_noMT.csv
│       └── liver_DE_HCC_vs_iCCA_noMT.csv
├── scripts/
│   └── ICB_reproduction_analysis.R
└── README.md
```
## Requirements

- R >= 4.0
- Seurat >= 5.0
- DoubletFinder
- dplyr
- ggplot2
- patchwork

## Usage
```r
# Run the complete analysis
source("scripts/ICB_reproduction_analysis.R")
```

## Key Findings

### Differential Expression Results
- **Melanoma**: 9,020 DE genes (post vs pre treatment)
- **BCC**: 5,210 DE genes (post vs pre treatment)
- **Liver**: 6,747 DE genes (HCC vs iCCA)

### Top Upregulated Genes Post-Treatment
**Melanoma**: PFN1, SERPINF1, EIF4A1, MAGEA4, IGF2BP1
**BCC**: RPL7A, RPS2, RARRES2, RAC3, MGST1
## Author
Trivedi Dhai
MS Bioinformatics, Northeastern University
December 2025

## References
Gondal et al. "Integrated cancer cell-specific single-cell RNA-seq datasets"

