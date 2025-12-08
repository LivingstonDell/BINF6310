
# ============================================
# 03_individual_visualizations.R
# ICB Reproduction Project - Individual Dataset Figures
# Author: Trivedi Dhai
# Date: November 2025
# ============================================

# DESCRIPTION:
# This script creates UMAP and QC visualizations for each 
# individual dataset before merging

# Run previous scripts first
source("/scratch/trivedi.dhai/ICB_reproduction/scripts/00_setup.R")
source("/scratch/trivedi.dhai/ICB_reproduction/scripts/01_data_loading.R")

# ============================================
# MELANOMA VISUALIZATIONS
# ============================================
print("Creating Melanoma figures...")
# UMAP by pre/post
p1 <- DimPlot(melanoma, reduction = "umap", group.by = "pre_post",
              cols = treatment_colors) +
  ggtitle("Melanoma: Pre vs Post ICB Treatment")
ggsave(paste0(paths$figures, "melanoma_umap_prepost.png"), plot = p1, width = 8, height = 6, dpi = 300)

# UMAP by patient
p2 <- DimPlot(melanoma, reduction = "umap", group.by = "patient") +
  ggtitle("Melanoma: By Patient")
ggsave(paste0(paths$figures, "melanoma_umap_patient.png"), plot = p2, width = 10, height = 6, dpi = 300)

# UMAP by study
p3 <- DimPlot(melanoma, reduction = "umap", group.by = "Study_name") +
  ggtitle("Melanoma: By Study")
ggsave(paste0(paths$figures, "melanoma_umap_study.png"), plot = p3, width = 8, height = 6, dpi = 300)

# QC violin plot
p4 <- VlnPlot(melanoma, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0)
ggsave(paste0(paths$figures, "melanoma_qc_violin.png"), plot = p4, width = 10, height = 5, dpi = 300)
# ============================================
# BCC VISUALIZATIONS
# ============================================
print("Creating BCC figures...")

# UMAP by pre/post
p5 <- DimPlot(bcc, reduction = "umap", group.by = "pre_post",
              cols = treatment_colors) +
  ggtitle("BCC: Pre vs Post ICB Treatment")
ggsave(paste0(paths$figures, "bcc_umap_prepost.png"), plot = p5, width = 8, height = 6, dpi = 300)

# UMAP by sample
p6 <- DimPlot(bcc, reduction = "umap", group.by = "sample_id") +
  ggtitle("BCC: By Patient")
ggsave(paste0(paths$figures, "bcc_umap_sample.png"), plot = p6, width = 8, height = 6, dpi = 300)

# QC violin plot
p7 <- VlnPlot(bcc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0)
ggsave(paste0(paths$figures, "bcc_qc_violin.png"), plot = p7, width = 10, height = 5, dpi = 300)
# ============================================
# LIVER VISUALIZATIONS
# ============================================
print("Creating Liver figures...")

# UMAP by cancer type
p8 <- DimPlot(liver, reduction = "umap", group.by = "Cancer_type",
              cols = c("HCC" = "#E41A1C", "iCCA" = "#377EB8")) +
  ggtitle("Liver: By Cancer Type")
ggsave(paste0(paths$figures, "liver_umap_cancertype.png"), plot = p8, width = 8, height = 6, dpi = 300)

# UMAP by pre/post
p9 <- DimPlot(liver, reduction = "umap", group.by = "pre_post",
              cols = treatment_colors) +
  ggtitle("Liver: Pre vs Post ICB Treatment")
ggsave(paste0(paths$figures, "liver_umap_prepost.png"), plot = p9, width = 8, height = 6, dpi = 300)

# UMAP by patient
p10 <- DimPlot(liver, reduction = "umap", group.by = "patient") +
  ggtitle("Liver: By Patient")
ggsave(paste0(paths$figures, "liver_umap_patient.png"), plot = p10, width = 10, height = 6, dpi = 300)

# QC violin plot
p11 <- VlnPlot(liver, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0)
ggsave(paste0(paths$figures, "liver_qc_violin.png"), plot = p11, width = 10, height = 5, dpi = 300)

# ============================================
# BREAST VISUALIZATIONS
# ============================================
print("Creating Breast figures...")

# UMAP by cancer subtype
p12 <- DimPlot(breast, reduction = "umap", group.by = "Cancer_type",
               cols = c("ER+" = "#E41A1C", "HER2+" = "#377EB8", "TNBC" = "#4DAF4A")) +
  ggtitle("Breast: By Cancer Subtype")
ggsave(paste0(paths$figures, "breast_umap_cancertype.png"), plot = p12, width = 8, height = 6, dpi = 300)

# UMAP by pre/post
p13 <- DimPlot(breast, reduction = "umap", group.by = "pre_post",
               cols = treatment_colors) +
  ggtitle("Breast: Pre vs Post ICB Treatment")
ggsave(paste0(paths$figures, "breast_umap_prepost.png"), plot = p13, width = 8, height = 6, dpi = 300)

# UMAP by outcome
p14 <- DimPlot(breast, reduction = "umap", group.by = "outcome",
               cols = c("E" = "#4DAF4A", "NE" = "#FF7F00", "n/a" = "#999999")) +
  ggtitle("Breast: By Outcome")
ggsave(paste0(paths$figures, "breast_umap_outcome.png"), plot = p14, width = 8, height = 6, dpi = 300)

# QC violin plot
p15 <- VlnPlot(breast, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0)
ggsave(paste0(paths$figures, "breast_qc_violin.png"), plot = p15, width = 10, height = 5, dpi = 300)

print("All individual dataset figures saved!")

