
# ============================================
# 05_combined_visualizations.R
# ICB Reproduction Project - Combined Dataset Figures
# Author: Trivedi Dhai
# Date: November 2025
# ============================================

# DESCRIPTION:
# This script creates visualizations for the merged dataset
# Reproduces figures similar to Gondal et al. Figure 3

# Run setup
source("/scratch/trivedi.dhai/ICB_reproduction/scripts/00_setup.R")

# Load processed combined data
combined <- readRDS(paste0(paths$processed_data, "combined_final.RDS"))
print(paste("Loaded combined object:", ncol(combined), "cells"))
# ============================================
# UMAP VISUALIZATIONS
# ============================================
print("Creating UMAP visualizations...")

# UMAP by dataset (Figure 3A style)
p1 <- DimPlot(combined, reduction = "umap", group.by = "dataset",
              cols = dataset_colors) +
  ggtitle("Combined ICB Datasets: Cancer Types")
ggsave(paste0(paths$figures, "combined_umap_dataset.png"), plot = p1, width = 10, height = 8, dpi = 300)

# UMAP by pre/post treatment
p2 <- DimPlot(combined, reduction = "umap", group.by = "pre_post",
              cols = treatment_colors) +
  ggtitle("Combined ICB Datasets: Pre vs Post Treatment")
ggsave(paste0(paths$figures, "combined_umap_prepost.png"), plot = p2, width = 10, height = 8, dpi = 300)

# UMAP split by dataset
p3 <- DimPlot(combined, reduction = "umap", group.by = "pre_post",
              split.by = "dataset", cols = treatment_colors) +
  ggtitle("Pre vs Post Treatment by Cancer Type")
ggsave(paste0(paths$figures, "combined_umap_split_dataset.png"), plot = p3, width = 16, height = 5, dpi = 300)

# ============================================
# PIE CHARTS (Figure 3C-E style)
# ============================================
print("Creating pie charts...")

# Dataset distribution
dataset_counts <- as.data.frame(table(combined$dataset))
colnames(dataset_counts) <- c("Dataset", "Cells")
dataset_counts$Percentage <- round(dataset_counts$Cells / sum(dataset_counts$Cells) * 100, 1)

p_pie1 <- ggplot(dataset_counts, aes(x = "", y = Cells, fill = Dataset)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  theme_void() +
  ggtitle("Cell Distribution Across Datasets") +
  scale_fill_manual(values = dataset_colors) +
  geom_text(aes(label = paste0(Percentage, "%")), 
            position = position_stack(vjust = 0.5), color = "white", size = 4)
ggsave(paste0(paths$figures, "combined_pie_dataset.png"), plot = p_pie1, width = 8, height = 6, dpi = 300)

# Pre/Post distribution
prepost_counts <- as.data.frame(table(combined$pre_post))
colnames(prepost_counts) <- c("Treatment", "Cells")
prepost_counts$Percentage <- round(prepost_counts$Cells / sum(prepost_counts$Cells) * 100, 1)

p_pie2 <- ggplot(prepost_counts, aes(x = "", y = Cells, fill = Treatment)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  theme_void() +
  ggtitle("Pre vs Post Treatment Distribution") +
  scale_fill_manual(values = treatment_colors) +
  geom_text(aes(label = paste0(Percentage, "%")), 
            position = position_stack(vjust = 0.5), color = "white", size = 5)
ggsave(paste0(paths$figures, "combined_pie_prepost.png"), plot = p_pie2, width = 8, height = 6, dpi = 300)

# ============================================
# BAR PLOTS
# ============================================
print("Creating bar plots...")

# Cell counts by dataset and treatment
count_data <- as.data.frame(table(combined$dataset, combined$pre_post))
colnames(count_data) <- c("Dataset", "Treatment", "Cells")

p_bar <- ggplot(count_data, aes(x = Dataset, y = Cells, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = treatment_colors) +
  theme_minimal() +
  ggtitle("Cell Counts by Dataset and Treatment Status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(paths$figures, "combined_barplot_treatment.png"), plot = p_bar, width = 10, height = 6, dpi = 300)

# ============================================
# QC COMPARISON
# ============================================
print("Creating QC comparison...")

p_qc <- VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA"), 
                group.by = "dataset", pt.size = 0, ncol = 2)
ggsave(paste0(paths$figures, "combined_qc_comparison.png"), plot = p_qc, width = 12, height = 5, dpi = 300)

# ============================================
# MARKER GENE EXPRESSION
# ============================================
print("Creating marker plots...")

# Epithelial markers
FeaturePlot(combined, features = c("EPCAM", "KRT18"), cols = c("lightgrey", "red"), ncol = 2)
ggsave(paste0(paths$figures, "combined_markers_epithelial.png"), width = 12, height = 5, dpi = 300)

VlnPlot(combined, features = c("EPCAM", "KRT18"), group.by = "dataset", pt.size = 0, ncol = 2)
ggsave(paste0(paths$figures, "combined_markers_violin.png"), width = 12, height = 5, dpi = 300)

print("All combined visualizations saved!")

