
# ============================================
# 07_de_visualizations.R
# ICB Reproduction Project - DE Visualization
# Author: Trivedi Dhai
# Date: November 2025
# ============================================

# DESCRIPTION:
# This script creates volcano plots, heatmaps, and other
# visualizations for differential expression results

# Run setup
source("/scratch/trivedi.dhai/ICB_reproduction/scripts/00_setup.R")

# Load data
combined <- readRDS(paste0(paths$processed_data, "combined_final.RDS"))
de_genes <- read.csv(paste0(paths$tables, "DE_post_vs_pre.csv"), row.names = 1)
de_melanoma <- read.csv(paste0(paths$tables, "DE_melanoma.csv"), row.names = 1)
de_bcc <- read.csv(paste0(paths$tables, "DE_bcc.csv"), row.names = 1)
de_liver <- read.csv(paste0(paths$tables, "DE_liver.csv"), row.names = 1)
de_breast <- read.csv(paste0(paths$tables, "DE_breast.csv"), row.names = 1)
# ============================================
# VOLCANO PLOT FUNCTION
# ============================================
create_volcano <- function(de_data, title) {
  de_data$gene <- rownames(de_data)
  de_data$significance <- "Not Significant"
  de_data$significance[de_data$avg_log2FC > 0.5 & de_data$p_val_adj < 0.05] <- "Up in Post"
  de_data$significance[de_data$avg_log2FC < -0.5 & de_data$p_val_adj < 0.05] <- "Down in Post"
  
  ggplot(de_data, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("Up in Post" = "#E41A1C", 
                                   "Down in Post" = "#377EB8", 
                                   "Not Significant" = "grey")) +
    theme_minimal() +
    ggtitle(title) +
    xlab("Log2 Fold Change") +
    ylab("-Log10 Adjusted P-value") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40")
}
# ============================================
# CREATE VOLCANO PLOTS
# ============================================
print("Creating volcano plots...")

p1 <- create_volcano(de_genes, "Combined: Post vs Pre ICB")
ggsave(paste0(paths$figures, "volcano_post_vs_pre.png"), plot = p1, width = 10, height = 8, dpi = 300)

p2 <- create_volcano(de_melanoma, "Melanoma: Post vs Pre ICB")
ggsave(paste0(paths$figures, "volcano_melanoma.png"), plot = p2, width = 8, height = 6, dpi = 300)

p3 <- create_volcano(de_bcc, "BCC: Post vs Pre ICB")
ggsave(paste0(paths$figures, "volcano_bcc.png"), plot = p3, width = 8, height = 6, dpi = 300)

p4 <- create_volcano(de_liver, "Liver: Post vs Pre ICB")
ggsave(paste0(paths$figures, "volcano_liver.png"), plot = p4, width = 8, height = 6, dpi = 300)

p5 <- create_volcano(de_breast, "Breast: Post vs Pre ICB")
ggsave(paste0(paths$figures, "volcano_breast.png"), plot = p5, width = 8, height = 6, dpi = 300)

# ============================================
# HEATMAPS
# ============================================
print("Creating heatmaps...")

# Get top DE genes
top_up <- rownames(de_genes[de_genes$avg_log2FC > 0, ])[1:20]
top_down <- rownames(de_genes[de_genes$avg_log2FC < 0, ])[1:20]
top_genes <- c(top_up, top_down)

# Scale genes for heatmap
combined <- ScaleData(combined, features = top_genes)

# Heatmap by treatment
DoHeatmap(combined, features = top_genes, group.by = "pre_post", size = 3) +
  ggtitle("Top 40 DE Genes: Post vs Pre ICB Treatment")
ggsave(paste0(paths$figures, "heatmap_top_DE_genes.png"), width = 12, height = 10, dpi = 300)

# Heatmap by dataset
DoHeatmap(combined, features = top_genes, group.by = "dataset", size = 3) +
  ggtitle("Top 40 DE Genes by Cancer Type")
ggsave(paste0(paths$figures, "heatmap_by_dataset.png"), width = 14, height = 10, dpi = 300)

# ============================================
# TOP GENE VISUALIZATION
# ============================================
print("Creating top gene visualizations...")

# Feature plot of top upregulated gene
FeaturePlot(combined, features = top_up[1], cols = c("lightgrey", "red")) +
  ggtitle(paste0("Top Upregulated Gene: ", top_up[1]))
ggsave(paste0(paths$figures, "top_upregulated_gene.png"), width = 8, height = 6, dpi = 300)

# Violin plot of top DE genes
VlnPlot(combined, features = c(top_up[1], top_down[1]), group.by = "pre_post",
        pt.size = 0, cols = treatment_colors)
ggsave(paste0(paths$figures, "top_DE_genes_violin.png"), width = 10, height = 5, dpi = 300)

# Save top genes table
top_genes_df <- rbind(
  data.frame(Direction = "Up_in_Post", gene = top_up, 
             avg_log2FC = de_genes[top_up, "avg_log2FC"],
             p_val_adj = de_genes[top_up, "p_val_adj"]),
  data.frame(Direction = "Down_in_Post", gene = top_down,
             avg_log2FC = de_genes[top_down, "avg_log2FC"],
             p_val_adj = de_genes[top_down, "p_val_adj"])
)
write.csv(top_genes_df, paste0(paths$tables, "top_DE_genes.csv"), row.names = FALSE)

print("All DE visualizations saved!")

