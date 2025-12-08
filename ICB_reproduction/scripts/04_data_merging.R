
# ============================================
# 04_data_merging.R
# ICB Reproduction Project - Data Merging
# Author: Trivedi Dhai
# Date: November 2025
# ============================================

# DESCRIPTION:
# This script merges all 4 datasets into a single Seurat object
# and performs normalization, scaling, PCA, and UMAP

# Run previous scripts first
source("/scratch/trivedi.dhai/ICB_reproduction/scripts/00_setup.R")
source("/scratch/trivedi.dhai/ICB_reproduction/scripts/01_data_loading.R")
# ============================================
# ADD DATASET LABELS
# ============================================
print("Adding dataset labels...")
melanoma$dataset <- "Melanoma"
bcc$dataset <- "BCC"
liver$dataset <- "Liver"
breast$dataset <- "Breast"

# ============================================
# MERGE DATASETS
# ============================================
print("Merging datasets...")
combined <- merge(melanoma, 
                  y = c(bcc, liver, breast),
                  add.cell.ids = c("Mel", "BCC", "Liv", "Br"))

print(paste("Combined object created:", ncol(combined), "cells"))
print("Distribution by dataset:")
print(table(combined$dataset))

# ============================================
# NORMALIZE DATA
# ============================================
print("Normalizing data...")
combined <- NormalizeData(combined)

# ============================================
# FIND VARIABLE FEATURES
# ============================================
print("Finding variable features...")
combined <- FindVariableFeatures(combined, 
                                  selection.method = "vst", 
                                  nfeatures = 2000)

# ============================================
# SCALE DATA
# ============================================
print("Scaling data...")
combined <- ScaleData(combined)
# ============================================
# RUN PCA
# ============================================
print("Running PCA...")
combined <- RunPCA(combined, npcs = 30)

# Elbow plot to determine dimensions
ElbowPlot(combined, ndims = 30)
ggsave(paste0(paths$figures, "combined_elbow_plot.png"), width = 8, height = 5, dpi = 300)

# ============================================
# RUN UMAP
# ============================================
print("Running UMAP...")
combined <- RunUMAP(combined, dims = 1:20)
# ============================================
# SAVE PROCESSED OBJECT
# ============================================
print("Saving combined object...")
saveRDS(combined, paste0(paths$processed_data, "combined_final.RDS"))

print("Data merging complete!")
print(combined)

