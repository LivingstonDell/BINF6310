
# ============================================
# 06_differential_expression.R
# ICB Reproduction Project - Differential Expression Analysis
# Author: Trivedi Dhai
# Date: November 2025
# ============================================

# DESCRIPTION:
# This script performs differential expression analysis
# comparing Post vs Pre ICB treatment cells

# Run setup
source("/scratch/trivedi.dhai/ICB_reproduction/scripts/00_setup.R")

# Load processed combined data
combined <- readRDS(paste0(paths$processed_data, "combined_final.RDS"))
print(paste("Loaded combined object:", ncol(combined), "cells"))
# ============================================
# COMBINED DE ANALYSIS (Post vs Pre)
# ============================================
print("Running combined DE analysis...")
Idents(combined) <- "pre_post"

de_genes <- FindMarkers(combined, 
                        ident.1 = "Post", 
                        ident.2 = "Pre",
                        max.cells.per.ident = 5000,
                        logfc.threshold = 0.25)

de_genes$gene <- rownames(de_genes)
print(paste("Combined DE genes found:", nrow(de_genes)))

# Save results
write.csv(de_genes, paste0(paths$tables, "DE_post_vs_pre.csv"))
# ============================================
# DATASET-SPECIFIC DE ANALYSIS
# ============================================

# Melanoma
print("Running Melanoma DE analysis...")
melanoma_cells <- subset(combined, dataset == "Melanoma")
Idents(melanoma_cells) <- "pre_post"
de_melanoma <- FindMarkers(melanoma_cells, ident.1 = "Post", ident.2 = "Pre", logfc.threshold = 0.25)
de_melanoma$gene <- rownames(de_melanoma)
write.csv(de_melanoma, paste0(paths$tables, "DE_melanoma.csv"))
print(paste("Melanoma DE genes:", nrow(de_melanoma)))

# BCC
print("Running BCC DE analysis...")
bcc_cells <- subset(combined, dataset == "BCC")
Idents(bcc_cells) <- "pre_post"
de_bcc <- FindMarkers(bcc_cells, ident.1 = "Post", ident.2 = "Pre", logfc.threshold = 0.25)
de_bcc$gene <- rownames(de_bcc)
write.csv(de_bcc, paste0(paths$tables, "DE_bcc.csv"))
print(paste("BCC DE genes:", nrow(de_bcc)))
# Liver
print("Running Liver DE analysis...")
liver_cells <- subset(combined, dataset == "Liver")
Idents(liver_cells) <- "pre_post"
de_liver <- FindMarkers(liver_cells, ident.1 = "Post", ident.2 = "Pre", logfc.threshold = 0.25)
de_liver$gene <- rownames(de_liver)
write.csv(de_liver, paste0(paths$tables, "DE_liver.csv"))
print(paste("Liver DE genes:", nrow(de_liver)))

# Breast
print("Running Breast DE analysis...")
breast_cells <- subset(combined, dataset == "Breast")
Idents(breast_cells) <- "pre_post"
de_breast <- FindMarkers(breast_cells, ident.1 = "Post", ident.2 = "Pre", 
                         max.cells.per.ident = 5000, logfc.threshold = 0.25)
de_breast$gene <- rownames(de_breast)
write.csv(de_breast, paste0(paths$tables, "DE_breast.csv"))
print(paste("Breast DE genes:", nrow(de_breast)))
# ============================================
# CREATE DE SUMMARY TABLE
# ============================================
de_summary <- data.frame(
  Dataset = c("Melanoma", "BCC", "Liver", "Breast", "Combined"),
  DE_Genes = c(nrow(de_melanoma), nrow(de_bcc), nrow(de_liver), nrow(de_breast), nrow(de_genes)),
  Total_Cells = c(ncol(melanoma_cells), ncol(bcc_cells), ncol(liver_cells), ncol(breast_cells), ncol(combined))
)

print("=== DE SUMMARY ===")
print(de_summary)
write.csv(de_summary, paste0(paths$tables, "DE_summary.csv"), row.names = FALSE)

print("Differential expression analysis complete!")

