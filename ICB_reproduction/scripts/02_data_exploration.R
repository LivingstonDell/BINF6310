 # ============================================
# 02_data_exploration.R
# ICB Reproduction Project - Data Exploration
# Author: Trivedi Dhai
# Date: November 2025
# ============================================

# DESCRIPTION:
# This script explores metadata and creates summary statistics
# for each individual dataset before merging

# Run previous scripts first
source("/scratch/trivedi.dhai/ICB_reproduction/scripts/00_setup.R")
source("/scratch/trivedi.dhai/ICB_reproduction/scripts/01_data_loading.R")

# ============================================
# EXPLORE MELANOMA DATASET
# ============================================
print("=== MELANOMA DATASET ===")
print(paste("Cells:", ncol(melanoma)))
print(paste("Genes:", nrow(melanoma)))
print("Metadata columns:")
print(colnames(melanoma@meta.data))
print("Pre/Post distribution:")
print(table(melanoma@meta.data$pre_post))
print("Patients:")
print(table(melanoma@meta.data$patient))

# ============================================
# EXPLORE BCC DATASET
# ============================================
print("=== BCC DATASET ===")
print(paste("Cells:", ncol(bcc)))
print(paste("Genes:", nrow(bcc)))
print("Pre/Post distribution:")
print(table(bcc@meta.data$pre_post))
print("Response:")
print(table(bcc@meta.data$Response))
print("Samples:")
print(table(bcc@meta.data$sample_id))

# ============================================
# EXPLORE LIVER DATASET
# ============================================
print("=== LIVER DATASET ===")
print(paste("Cells:", ncol(liver)))
print(paste("Genes:", nrow(liver)))
print("Cancer types (HCC vs iCCA):")
print(table(liver@meta.data$Cancer_type))
print("Pre/Post distribution:")
print(table(liver@meta.data$pre_post))
print("Patients:")
print(table(liver@meta.data$patient))

# ============================================
# EXPLORE BREAST DATASET
# ============================================
print("=== BREAST DATASET ===")
print(paste("Cells:", ncol(breast)))
print(paste("Genes:", nrow(breast)))
print("Cancer subtypes:")
print(table(breast@meta.data$Cancer_type))
print("Pre/Post distribution:")
print(table(breast@meta.data$pre_post))
print("Outcome:")
print(table(breast@meta.data$outcome))
# ============================================
# CREATE SUMMARY TABLE
# ============================================
summary_df <- data.frame(
  Dataset = c("Melanoma", "BCC", "Liver", "Breast"),
  Study = c("Jerby-Arnon/Tirosh", "Yost", "Ma", "Bassez"),
  Cancer_Type = c("Skin", "Skin", "Liver (HCC+iCCA)", "Breast (ER+, HER2+, TNBC)"),
  Total_Cells = c(ncol(melanoma), ncol(bcc), ncol(liver), ncol(breast)),
  Pre_Treatment = c(sum(melanoma@meta.data$pre_post == "Pre"),
                    sum(bcc@meta.data$pre_post == "Pre"),
                    sum(liver@meta.data$pre_post == "Pre"),
                    sum(breast@meta.data$pre_post == "Pre")),
  Post_Treatment = c(sum(melanoma@meta.data$pre_post == "Post"),
                     sum(bcc@meta.data$pre_post == "Post"),
                     sum(liver@meta.data$pre_post == "Post"),
                     sum(breast@meta.data$pre_post == "Post"))
)

print("=== SUMMARY TABLE ===")
print(summary_df)

# Save summary
write.csv(summary_df, paste0(paths$tables, "dataset_summary.csv"), row.names = FALSE)
print("Summary table saved!")

