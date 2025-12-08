# ICB Raw Data Reproduction Project

## Integrated Cancer Cell-Specific Single-Cell RNA-seq Analysis from Raw GEO Data

**Authors:** Dhaivat Trivedi, Sneha Kini, Viswa Varenya, Trang Do  
**Course:** BINF 6310: Final Project  
**Institution:** Northeastern University  
**Date:** December 2025

---

## 📖 Project Overview

This project reproduces the single-cell RNA-seq analysis from Gondal et al. (2025) using **raw data downloaded directly from GEO**, rather than the authors' pre-processed RDS objects. This approach provides a more rigorous test of reproducibility by implementing the complete analysis pipeline from scratch.

**Objective:** Validate the reproducibility of the original study's methodology by processing raw count matrices through the full scRNA-seq pipeline including QC, DoubletFinder, normalization, and differential expression analysis.

---

## 📄 Original Publication

**Title:** Integrated cancer cell-specific single-cell RNA-seq datasets of immune checkpoint blockade-treated patients

**Authors:** Gondal, M.N., Cieslik, M. & Chinnaiyan, A.M.

**Journal:** Scientific Data (Nature)

**Year:** 2025

**DOI:** [10.1038/s41597-025-04381-6](https://doi.org/10.1038/s41597-025-04381-6)

**GitHub:** [github.com/LivingstonDell/BINF6310](https://github.com/LivingstonDell/BINF6310)

---

## 🗂️ Datasets

| Dataset | GEO Accession | Cancer Type | Our Cells | Paper Cells | Match |
|---------|---------------|-------------|-----------|-------------|-------|
| Melanoma | GSE115978 | Skin melanoma | 2,018 | 1,945 | 96% |
| BCC | GSE123813 | Basal cell carcinoma | 3,452 | 3,500 | 99% |
| Liver | GSE125449 | HCC/iCCA | 1,577 | 1,992 | 79% |

---

## 🔬 Analysis Pipeline

| Step | Description | Tools Used |
|------|-------------|------------|
| 1. Data Download | Raw counts from GEO | wget, curl |
| 2. Seurat Object Creation | Load count matrices | CreateSeuratObject() |
| 3. Cell Type Filtering | Keep malignant/tumor cells only | Metadata filtering |
| 4. DoubletFinder | Remove doublets (10X data only) | DoubletFinder package |
| 5. Standard Processing | Normalize, Scale, PCA, UMAP, Clustering | Seurat pipeline |
| 6. Differential Expression | Pre vs Post treatment comparison | FindMarkers() |
| 7. Visualization | UMAP and Volcano plots | ggplot2 |

---

## 📁 Directory Structure

```
ICB_raw_reproduction/
├── README.md                              # This file
├── data/
│   ├── raw/                               # Raw data from GEO
│   │   ├── melanoma/
│   │   │   ├── GSE115978_counts.csv.gz
│   │   │   └── GSE115978_cell.annotations.csv.gz
│   │   ├── bcc/
│   │   │   └── GSE123813_*.mtx.gz
│   │   └── liver/
│   │       └── GSE125449_*.mtx.gz
│   └── processed/                         # Processed Seurat objects
│       ├── melanoma_final_processed.RDS
│       ├── bcc_final_processed.RDS
│       └── liver_final_processed.RDS
├── scripts/
│   └── ICB_reproduction_analysis.R        # Main analysis script
├── results/
│   ├── figures/                           # Output plots
│   │   ├── melanoma_umap.png
│   │   ├── bcc_umap.png
│   │   ├── liver_umap.png
│   │   ├── volcano_melanoma.png
│   │   ├── volcano_bcc.png
│   │   └── volcano_liver.png
│   └── tables/                            # Output tables
│       ├── melanoma_DE_post_vs_pre.csv
│       ├── bcc_DE_post_vs_pre_noMT.csv
│       └── liver_DE_HCC_vs_iCCA_noMT.csv
└── docs/                                  # Documentation
```

---

## 📊 Results Summary

### Figures Generated

| Figure | Description |
|--------|-------------|
| `melanoma_umap.png` | Melanoma UMAP colored by treatment status |
| `bcc_umap.png` | BCC UMAP colored by treatment status |
| `liver_umap.png` | Liver UMAP colored by cancer type (HCC vs iCCA) |
| `volcano_melanoma.png` | Volcano plot - Melanoma DE genes |
| `volcano_bcc.png` | Volcano plot - BCC DE genes |
| `volcano_liver.png` | Volcano plot - Liver DE genes |

### Tables Generated

| Table | Description |
|-------|-------------|
| `melanoma_DE_post_vs_pre.csv` | Melanoma DE results (9,020 genes) |
| `bcc_DE_post_vs_pre_noMT.csv` | BCC DE results (5,210 genes) |
| `liver_DE_HCC_vs_iCCA_noMT.csv` | Liver DE results (6,747 genes) |

---

## 🔬 Key Findings

### Differential Expression Summary

| Dataset | DE Genes (Post vs Pre) | Our Cells | Paper Cells |
|---------|------------------------|-----------|-------------|
| Melanoma | 9,020 | 2,018 | 1,945 |
| BCC | 5,210 | 3,452 | 3,500 |
| Liver | 6,747 | 1,577 | 1,992 |

### Top Upregulated Genes Post-Treatment

**Melanoma:**
| Gene | Function |
|------|----------|
| PFN1 | Cytoskeleton remodeling, cell migration |
| SERPINF1 | Anti-angiogenic factor |
| MAGEA4 | Melanoma tumor antigen |
| IGF2BP1 | Promotes tumor progression |

**BCC:**
| Gene | Function |
|------|----------|
| RPL7A, RPS2, RPL5 | Ribosomal proteins (protein synthesis) |
| RARRES2 | Immune cell recruitment |
| RAC3, MGST1 | Cell migration, oxidative stress |

### Biological Interpretation

Cancer cells surviving ICB treatment show:
- **Increased protein synthesis** — stress response mechanism
- **Enhanced cell motility** — potential immune evasion strategy
- **Translation machinery upregulation** — consistent across skin cancers

---

## 🚀 How to Run

### Prerequisites
- R version 4.0+
- Seurat >= 5.0
- DoubletFinder
- dplyr
- ggplot2
- patchwork

### Setup Environment (HPC)
```bash
# Get interactive session
srun --pty --mem=16G --time=03:00:00 bash

# Activate conda environment
conda activate icb_analysis
module unload R

# Navigate to project
cd /scratch/$USER/ICB_raw_reproduction

# Start R
R
```

### Run Complete Analysis
```r
# Run the complete analysis script
source("scripts/ICB_reproduction_analysis.R")
```

### Load Pre-processed Data
```r
# Load pre-processed objects
library(Seurat)
library(dplyr)
library(ggplot2)

processed_path <- "data/processed/"
melanoma <- readRDS(paste0(processed_path, "melanoma_final_processed.RDS"))
bcc <- readRDS(paste0(processed_path, "bcc_final_processed.RDS"))
liver <- readRDS(paste0(processed_path, "liver_final_processed.RDS"))
```

---

## 📦 Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| R | 4.4.3 | Base language |
| Seurat | 5.0+ | Single-cell analysis |
| DoubletFinder | latest | Doublet detection (10X data) |
| dplyr | 1.1+ | Data manipulation |
| ggplot2 | 3.4+ | Visualization |
| patchwork | 1.1+ | Plot arrangement |

### Installation
```r
install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork"))
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
```

---

## 🔗 Data Sources

| Source | Accession | URL |
|--------|-----------|-----|
| GEO - Melanoma | GSE115978 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115978 |
| GEO - BCC | GSE123813 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123813 |
| GEO - Liver | GSE125449 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125449 |
| GitHub (Our Code) | - | https://github.com/LivingstonDell/BINF6310 |

---

## 📈 Reproducibility Notes

### What Was Successfully Reproduced
- ✅ Raw data download and processing from GEO
- ✅ DoubletFinder pipeline for 10X data (BCC, Liver)
- ✅ Cell type filtering to malignant cells only
- ✅ UMAP visualization and clustering
- ✅ Differential expression analysis
- ✅ Volcano plot generation

### Cell Count Comparison
| Dataset | Our Count | Paper Count | Difference | Notes |
|---------|-----------|-------------|------------|-------|
| Melanoma | 2,018 | 1,945 | +73 (96%) | Excellent match |
| BCC | 3,452 | 3,500 | -48 (99%) | Excellent match |
| Liver | 1,577 | 1,992 | -415 (79%) | DoubletFinder parameters |

### Critical Insights Identified

**Issue 1: DoubletFinder Parameter Sensitivity**
- Default 0.8% doublet rate may not fit all datasets
- Liver showed larger discrepancy due to parameter choices

**Issue 2: Metadata Dependency**
- Liver required external annotation files to identify malignant cells
- Without metadata: 9,150 cells → With filtering: 1,577 cells

**Issue 3: MT Gene Filtering**
- Initial Liver DE dominated by mitochondrial genes
- Paper doesn't explicitly state MT gene removal criteria

### Recommendations for Future Reproducibility
- Use Harmony or scVI for batch correction instead of hSEGs
- Implement Scrublet + DoubletFinder consensus for doublet detection
- Document cell counts at each filtering step

---

## 📚 References

1. Gondal, M.N., Cieslik, M. & Chinnaiyan, A.M. Integrated cancer cell-specific single-cell RNA-seq datasets of immune checkpoint blockade-treated patients. *Sci Data* **12**, 139 (2025).

2. Jerby-Arnon, L. et al. A cancer cell program promotes T cell exclusion and resistance to checkpoint blockade. *Cell* **175**, 984–997 (2018).

3. Yost, K.E. et al. Clonal replacement of tumor-specific T cells following PD-1 blockade. *Nat Med* **25**, 1251–1259 (2019).

4. Ma, L. et al. Tumor cell biodiversity drives microenvironmental reprogramming in liver cancer. *Cancer Cell* **36**, 418–430 (2019).

5. McGinnis, C.S., Murrow, L.M. & Gartner, Z.J. DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors. *Cell Syst* **8**, 329–337 (2019).

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

*Last updated: December 2025*
