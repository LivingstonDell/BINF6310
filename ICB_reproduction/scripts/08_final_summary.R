 
# ============================================
# 08_final_summary.R
# ICB Reproduction Project - Final Summary
# Author: Trivedi Dhai
# Date: November 2025
# ============================================

# DESCRIPTION:
# This script generates final summary statistics and project logs

# Run setup
source("/scratch/trivedi.dhai/ICB_reproduction/scripts/00_setup.R")

# Load data
combined <- readRDS(paste0(paths$processed_data, "combined_final.RDS"))

# ============================================
# COMPREHENSIVE SUMMARY TABLE
# ============================================
summary_full <- data.frame(
  Dataset = c("Melanoma", "BCC", "Liver", "Breast", "TOTAL"),
  Study = c("Jerby-Arnon/Tirosh", "Yost", "Ma", "Bassez", "-"),
  Total_Cells = c(
    sum(combined$dataset == "Melanoma"),
    sum(combined$dataset == "BCC"),
    sum(combined$dataset == "Liver"),
    sum(combined$dataset == "Breast"),
    ncol(combined)
  ),
  Pre_Treatment = c(
    sum(combined$dataset == "Melanoma" & combined$pre_post == "Pre"),
    sum(combined$dataset == "BCC" & combined$pre_post == "Pre"),
    sum(combined$dataset == "Liver" & combined$pre_post == "Pre"),
    sum(combined$dataset == "Breast" & combined$pre_post == "Pre"),
    sum(combined$pre_post == "Pre")
  ),
  Post_Treatment = c(
    sum(combined$dataset == "Melanoma" & combined$pre_post == "Post"),
    sum(combined$dataset == "BCC" & combined$pre_post == "Post"),
    sum(combined$dataset == "Liver" & combined$pre_post == "Post"),
    sum(combined$dataset == "Breast" & combined$pre_post == "Post"),
    sum(combined$pre_post == "Post")
  )
)
print("=== COMPREHENSIVE SUMMARY ===")
print(summary_full)
write.csv(summary_full, paste0(paths$tables, "comprehensive_summary.csv"), row.names = FALSE)

# ============================================
# PROJECT LOG
# ============================================
project_log <- data.frame(
  Item = c("Total Cells Analyzed",
           "Total Genes",
           "Datasets",
           "Figures Generated",
           "Tables Created",
           "Analysis Date"),
  Value = c(ncol(combined),
            nrow(combined),
            4,
            length(list.files(paths$figures)),
            length(list.files(paths$tables)),
            as.character(Sys.Date()))
)

print("=== PROJECT LOG ===")
print(project_log)
write.csv(project_log, paste0(paths$tables, "project_log.csv"), row.names = FALSE)

# ============================================
# LIST ALL OUTPUTS
# ============================================
print("=== FIGURES GENERATED ===")
print(list.files(paths$figures))

print("=== TABLES GENERATED ===")
print(list.files(paths$tables))

print("=== SCRIPTS ===")
print(list.files(paths$scripts))

print("Project summary complete!")

