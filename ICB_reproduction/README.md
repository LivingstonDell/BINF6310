# ICB scRNA-seq Reproducibility Project

## Integrated Cancer Cell-Specific Single-Cell RNA-seq Analysis of ICB-Treated Patients

**Authors:** Dhaivat Trivedi, Sneha Kini, Viswa Varenya, Trang Do  
**Course:** BINF 6310: Final Project  
**Institution:** Northeastern University  
**Date:** November 2025

---

## 📖 Project Overview

This project reproduces the computational analysis from Gondal et al. (2025), which integrated cancer cell-specific single-cell RNA-seq datasets from immune checkpoint blockade (ICB) treated patients across multiple cancer types.

**Objective:** Validate the reproducibility of the original study's findings using the authors' processed datasets and standard scRNA-seq analysis pipelines.

---

## 📄 Original Publication

**Title:** Integrated cancer cell-specific single-cell RNA-seq datasets of immune checkpoint blockade-treated patients

**Authors:** Gondal, M.N., Cieslik, M. & Chinnaiyan, A.M.

**Journal:** Scientific Data (Nature)

**Year:** 2025

**DOI:** [10.1038/s41597-025-04381-6](https://doi.org/10.1038/s41597-025-04381-6)

**Data Repository:** [Zenodo - 10.5281/zenodo.10407125](https://doi.org/10.5281/zenodo.10407125)

---

## 🗂️ Datasets

| Dataset | Study | Cancer Type | Total Cells | Pre-Treatment | Post-Treatment |
|---------|-------|-------------|-------------|---------------|----------------|
| Melanoma | Jerby-Arnon/Tirosh | Skin | 1,945 | 1,143 | 802 |
| BCC | Yost et al. | Basal Cell Carcinoma | 3,500 | 2,521 | 979 |
| Liver | Ma et al. | HCC + iCCA | 1,992 | 275 | 1,717 |
| Breast | Bassez et al. | TNBC, HER2+, ER+ | 65,882 | 35,587 | 30,295 |
| **TOTAL** | - | - | **73,319** | **39,526** | **33,793** |

---

## 📁 Directory Structure

```
ICB_reproduction/
├── README.md                          # This file
├── data/
│   ├── raw/                           # Original RDS files from Zenodo
│   │   ├── jerby_arnon_melanoma.RDS
│   │   ├── yost_bcc.RDS
│   │   ├── ma_liver.RDS
│   │   └── bassez_breast.RDS
│   └── processed/                     # Processed Seurat objects
│       ├── melanoma_processed.RDS
│       ├── bcc_processed.RDS
│       ├── liver_processed.RDS
│       ├── breast_processed.RDS
│       └── combined_final.RDS
├── scripts/                           # R analysis scripts (00-08)
│   ├── 00_setup.R
│   ├── 01_data_loading.R
│   ├── 02_data_exploration.R
│   ├── 03_individual_visualizations.R
│   ├── 04_data_merging.R
│   ├── 05_combined_visualizations.R
│   ├── 06_differential_expression.R
│   ├── 07_de_visualizations.R
│   └── 08_final_summary.R
├── results/
│   ├── figures/                       # Output plots (34 PNG files)
│   └── tables/                        # Output tables (10 CSV files)
└── docs/                              # Documentation
```

---

## 📜 Scripts Description

| Script | Description | Key Outputs |
|--------|-------------|-------------|
| `00_setup.R` | Load libraries, define paths, set color palettes | Environment setup |
| `01_data_loading.R` | Load all 4 RDS datasets from Zenodo | Seurat objects |
| `02_data_exploration.R` | Explore metadata, create summary statistics | dataset_summary.csv |
| `03_individual_visualizations.R` | UMAPs and QC plots for each dataset | 15 individual figures |
| `04_data_merging.R` | Merge datasets, normalize, PCA, UMAP | combined_final.RDS |
| `05_combined_visualizations.R` | Combined UMAPs, pie charts, bar plots | 8 combined figures |
| `06_differential_expression.R` | DE analysis (Post vs Pre) for all datasets | DE tables (5 CSV) |
| `07_de_visualizations.R` | Volcano plots, heatmaps, top gene plots | 6 DE figures |
| `08_final_summary.R` | Generate final summary and project log | comprehensive_summary.csv |

---

## 📊 Results Summary

### Figures Generated (34 total)

#### Individual Dataset Figures
| Figure | Description |
|--------|-------------|
| `melanoma_umap_prepost.png` | Melanoma UMAP colored by treatment status |
| `melanoma_umap_patient.png` | Melanoma UMAP colored by patient |
| `melanoma_umap_study.png` | Melanoma UMAP colored by study origin |
| `melanoma_qc_violin.png` | Melanoma QC metrics violin plot |
| `bcc_umap_prepost.png` | BCC UMAP colored by treatment status |
| `bcc_umap_sample.png` | BCC UMAP colored by sample |
| `bcc_qc_violin.png` | BCC QC metrics violin plot |
| `liver_umap_cancertype.png` | Liver UMAP colored by cancer type (HCC vs iCCA) |
| `liver_umap_prepost.png` | Liver UMAP colored by treatment status |
| `liver_umap_patient.png` | Liver UMAP colored by patient |
| `liver_qc_violin.png` | Liver QC metrics violin plot |
| `breast_umap_cancertype.png` | Breast UMAP colored by subtype |
| `breast_umap_prepost.png` | Breast UMAP colored by treatment status |
| `breast_umap_outcome.png` | Breast UMAP colored by outcome |
| `breast_qc_violin.png` | Breast QC metrics violin plot |

#### Combined Dataset Figures
| Figure | Description |
|--------|-------------|
| `combined_elbow_plot.png` | PCA elbow plot for dimension selection |
| `combined_umap_dataset.png` | Combined UMAP colored by cancer type |
| `combined_umap_prepost.png` | Combined UMAP colored by treatment |
| `combined_umap_split_dataset.png` | Split UMAP by dataset |
| `combined_pie_dataset.png` | Pie chart of cell distribution |
| `combined_pie_prepost.png` | Pie chart of treatment distribution |
| `combined_barplot_treatment.png` | Bar plot of cells by treatment |
| `combined_qc_comparison.png` | QC comparison across datasets |
| `combined_markers_epithelial.png` | Epithelial marker expression |
| `combined_markers_violin.png` | Marker violin plots |

#### Differential Expression Figures
| Figure | Description |
|--------|-------------|
| `volcano_post_vs_pre.png` | Volcano plot - combined analysis |
| `volcano_melanoma.png` | Volcano plot - melanoma |
| `volcano_bcc.png` | Volcano plot - BCC |
| `volcano_liver.png` | Volcano plot - liver |
| `volcano_breast.png` | Volcano plot - breast |
| `heatmap_top_DE_genes.png` | Heatmap of top 40 DE genes |
| `heatmap_by_dataset.png` | Heatmap split by dataset |
| `top_upregulated_gene.png` | Top upregulated gene feature plot |
| `top_DE_genes_violin.png` | Violin plot of top DE genes |

### Tables Generated (10 total)

| Table | Description |
|-------|-------------|
| `dataset_summary.csv` | Basic dataset statistics |
| `comprehensive_summary.csv` | Full summary with pre/post counts |
| `DE_post_vs_pre.csv` | Combined DE results (7,767 genes) |
| `DE_melanoma.csv` | Melanoma-specific DE (13,649 genes) |
| `DE_bcc.csv` | BCC-specific DE (9,909 genes) |
| `DE_liver.csv` | Liver-specific DE (13,214 genes) |
| `DE_breast.csv` | Breast-specific DE (6,754 genes) |
| `DE_summary.csv` | Summary of DE gene counts |
| `top_DE_genes.csv` | Top 40 DE genes (20 up, 20 down) |
| `project_log.csv` | Project completion log |

---

## 🔬 Key Findings

### Differential Expression Summary

| Dataset | DE Genes (Post vs Pre) | Total Cells |
|---------|------------------------|-------------|
| Melanoma | 13,649 | 1,945 |
| BCC | 9,909 | 3,500 |
| Liver | 13,214 | 1,992 |
| Breast | 6,754 | 65,882 |
| **Combined** | **7,767** | **73,319** |

### Top Upregulated Genes (Post-ICB)
- IGF2 (log2FC: +1.65)
- SERPINA1 (log2FC: +1.78)
- CHGB (log2FC: +3.58)

### Top Downregulated Genes (Post-ICB)
- IGKV4-1 (log2FC: -2.61) - Immunoglobulin
- IGKV3-20 (log2FC: -1.67) - Immunoglobulin
- HLA-B (log2FC: -0.78) - Immune presentation

---

## 🚀 How to Run

### Prerequisites
- R version 4.3+
- Seurat package
- dplyr
- ggplot2

### Setup Environment (HPC)
```bash
# Get interactive session
srun --pty --mem=16G --time=02:00:00 bash

# Activate conda environment
conda activate icb_analysis
module unload R

# Start R
R
```

### Run Scripts Sequentially
```r
# Option 1: Run all scripts in order
source("scripts/00_setup.R")
source("scripts/01_data_loading.R")
source("scripts/02_data_exploration.R")
source("scripts/03_individual_visualizations.R")
source("scripts/04_data_merging.R")
source("scripts/05_combined_visualizations.R")
source("scripts/06_differential_expression.R")
source("scripts/07_de_visualizations.R")
source("scripts/08_final_summary.R")
```

### Load Pre-processed Data
```r
# Option 2: Load pre-processed combined object
library(Seurat)
library(dplyr)
library(ggplot2)

combined <- readRDS("data/processed/combined_final.RDS")
de_genes <- read.csv("results/tables/DE_post_vs_pre.csv", row.names = 1)
```

---

## 📦 Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| R | 4.4.3 | Base language |
| Seurat | 5.0+ | Single-cell analysis |
| dplyr | 1.1+ | Data manipulation |
| ggplot2 | 3.4+ | Visualization |
| SeuratObject | 5.0+ | Seurat data structures |

### Installation
```r
install.packages(c("Seurat", "dplyr", "ggplot2"))
```

---

## 🔗 Data Sources

| Source | URL |
|--------|-----|
| Zenodo (Processed RDS) | https://doi.org/10.5281/zenodo.10407125 |
| GitHub (Original Code) | https://github.com/MahnoorNGondal/scRNA-seq-ICB-cohorts |
| GEO - Melanoma | GSE115978 |
| GEO - BCC | GSE123813 |
| GEO - Liver | GSE125449 |

---

## 📈 Reproducibility Notes

### What Was Successfully Reproduced
- ✅ Data loading and integration from 4 cancer types
- ✅ Combined UMAP visualization showing cancer type separation
- ✅ Pre vs Post treatment differential expression analysis
- ✅ Cell count distributions matching paper's Figure 3

### Minor Differences Observed
- Exact cell counts may vary slightly due to software version differences
- Color schemes adapted for clarity
- Additional QC visualizations added

### Limitations
- Used authors' pre-processed RDS files (not raw count matrices)
- Did not reproduce all supplementary analyses
- DoubletFinder QC not applied (data pre-filtered)

---

## 📚 References

1. Gondal, M.N., Cieslik, M. & Chinnaiyan, A.M. Integrated cancer cell-specific single-cell RNA-seq datasets of immune checkpoint blockade-treated patients. *Sci Data* **12**, 139 (2025).

2. Jerby-Arnon, L. et al. A cancer cell program promotes T cell exclusion and resistance to checkpoint blockade. *Cell* **175**, 984–997 (2018).

3. Yost, K.E. et al. Clonal replacement of tumor-specific T cells following PD-1 blockade. *Nat Med* **25**, 1251–1259 (2019).

4. Ma, L. et al. Tumor cell biodiversity drives microenvironmental reprogramming in liver cancer. *Cancer Cell* **36**, 418–430 (2019).

5. Bassez, A. et al. A single-cell map of intratumoral changes during anti-PD1 treatment of patients with breast cancer. *Nat Med* **27**, 820–832 (2021).

---

## 📝 License

This project is the Final Project for BINF 6310 at Northeastern University.

---

## 👥 Authors

- Dhaivat Trivedi
- Sneha Kini
- Viswa Varenya
- Trang Do

**Course:** BINF 6310: Final Project  
**Institution:** Northeastern University

---

*Last updated: November 2025*
