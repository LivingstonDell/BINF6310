# ============================================
# 01_data_loading.R
# ICB Reproduction Project - Data Loading
# Author: Trivedi Dhai
# Date: November 2025
# ============================================

# DESCRIPTION:
# This script loads all 4 datasets from the ICB study
# Datasets: Melanoma, BCC, Liver (HCC+iCCA), Breast (TNBC, HER2+, ER+)
# Source: Gondal et al., Scientific Data (2025)

# Run 00_setup.R first
source("/scratch/trivedi.dhai/ICB_reproduction/scripts/00_setup.R")

# ============================================
# LOAD DATASETS
# ============================================

# Dataset 1: Melanoma (Jerby-Arnon & Tirosh)
# Cancer type: Skin melanoma
# Technology: SmartSeq2
# Cells: ~1,945
melanoma <- readRDS(paste0(paths$raw_data, "jerby_arnon_melanoma.RDS"))
print(paste("Melanoma loaded:", ncol(melanoma), "cells"))

# Dataset 2: BCC (Yost et al.)
# Cancer type: Basal Cell Carcinoma
# Technology: 10X
# Cells: ~3,500
bcc <- readRDS(paste0(paths$raw_data, "yost_bcc.RDS"))
print(paste("BCC loaded:", ncol(bcc), "cells"))

# Dataset 3: Liver (Ma et al.)
# Cancer type: HCC + iCCA
# Technology: 10X
# Cells: ~1,992
liver <- readRDS(paste0(paths$raw_data, "ma_liver.RDS"))
print(paste("Liver loaded:", ncol(liver), "cells"))

# Dataset 4: Breast (Bassez et al.)
# Cancer type: TNBC, HER2+, ER+
# Technology: 10X
# Cells: ~65,882
breast <- readRDS(paste0(paths$raw_data, "bassez_breast.RDS"))
print(paste("Breast loaded:", ncol(breast), "cells"))

# ============================================
# VERIFY DATASETS
# ============================================
print("=== Dataset Summary ===")
print(paste("Total cells:", ncol(melanoma) + ncol(bcc) + ncol(liver) + ncol(breast)))

# View structure of each dataset
print("Melanoma structure:")
print(melanoma)

print("BCC structure:")
print(bcc)

print("Liver structure:")
print(liver)

print("Breast structure:")
print(breast)
print(paste("Total cells:", ncol(melanoma) + ncol(bcc) + ncol(liver) + ncol(breast)))

# View structure of each dataset
print("Melanoma structure:")
print(melanoma)

print("BCC structure:")
print(bcc)

print("Liver structure:")
print(liver)

print("Breast structure:")
print(breast)
print(paste("Total cells:", ncol(melanoma) + ncol(bcc) + ncol(liver) + ncol(breast)))
# View structure of each dataset
print("Melanoma structure:")
print(melanoma)

print("BCC structure:")
print(bcc)

print("Liver structure:")
print(liver)

print("Breast structure:")
print(breast)

