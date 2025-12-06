
# ============================================
# ICB Raw Data Reproduction Analysis
# Author: Trivedi Dhai
# Date: December 2025
# ============================================

# This script reproduces the analysis from:
# "Integrated cancer cell-specific single-cell RNA-seq datasets"
# Using raw data from GEO

# ============================================
# 1. SETUP AND LIBRARIES
# ============================================

library(Seurat)
library(DoubletFinder)
library(dplyr)
library(ggplot2)
library(patchwork)

# Set paths
base_path <- "/scratch/trivedi.dhai/ICB_raw_reproduction/"
raw_path <- paste0(base_path, "data/raw/")
processed_path <- paste0(base_path, "data/processed/")
fig_path <- paste0(base_path, "results/figures/")
table_path <- paste0(base_path, "results/tables/")

# Create directories
dir.create(processed_path, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)
dir.create(table_path, recursive = TRUE, showWarnings = FALSE)
# ============================================
# 2. LOAD RAW DATA
# ============================================

# --- MELANOMA (GSE115978) ---
melanoma_counts <- read.table(gzfile(paste0(raw_path, "melanoma/GSE115978_counts.csv.gz")), 
                              header = TRUE, sep = ",", row.names = 1)
melanoma_meta <- read.table(gzfile(paste0(raw_path, "melanoma/GSE115978_cell.annotations.csv.gz")), 
                            header = TRUE, sep = ",", row.names = 1)
melanoma <- CreateSeuratObject(counts = melanoma_counts, meta.data = melanoma_meta, project = "Melanoma")

# --- BCC (GSE123813) ---
bcc_counts <- Read10X(data.dir = paste0(raw_path, "bcc/"))
bcc_meta <- read.table(gzfile(paste0(raw_path, "bcc/GSE123813_bcc_scRNA_counts.txt.gz")), 
                       header = TRUE, nrows = 1)
bcc <- CreateSeuratObject(counts = bcc_counts, project = "BCC")

# --- LIVER (GSE125449) ---
# Set 1
set1_barcodes <- read.table(gzfile(paste0(raw_path, "liver/GSE125449_Set1_barcodes.tsv.gz")), header = FALSE)
set1_genes <- read.table(gzfile(paste0(raw_path, "liver/GSE125449_Set1_genes.tsv.gz")), header = FALSE)
set1_matrix <- Matrix::readMM(gzfile(paste0(raw_path, "liver/GSE125449_Set1_matrix.mtx.gz")))
rownames(set1_matrix) <- set1_genes$V1
colnames(set1_matrix) <- set1_barcodes$V1

# Set 2
set2_barcodes <- read.table(gzfile(paste0(raw_path, "liver/GSE125449_Set2_barcodes.tsv.gz")), header = FALSE)
set2_genes <- read.table(gzfile(paste0(raw_path, "liver/GSE125449_Set2_genes.tsv.gz")), header = FALSE)
set2_matrix <- Matrix::readMM(gzfile(paste0(raw_path, "liver/GSE125449_Set2_matrix.mtx.gz")))
rownames(set2_matrix) <- set2_genes$V1
colnames(set2_matrix) <- set2_barcodes$V1

liver <- CreateSeuratObject(counts = cbind(set1_matrix, set2_matrix), project = "Liver")
# ============================================
# 3. FILTER TO CANCER CELLS ONLY
# ============================================

# Melanoma - filter to "Mal" (malignant) cells
melanoma_cancer <- subset(melanoma, cell.types == "Mal")

# BCC - filter to "Tumor" cells
bcc_cancer <- subset(bcc, cluster %in% grep("Tumor", unique(bcc$cluster), value = TRUE))

# Liver - filter using sample files (will be done after DoubletFinder)

# ============================================
# 4. DOUBLETFINDER (10X data only)
# ============================================

# --- BCC DoubletFinder ---
bcc_cancer <- NormalizeData(bcc_cancer)
bcc_cancer <- FindVariableFeatures(bcc_cancer, selection.method = "vst", nfeatures = 2000)
bcc_cancer <- ScaleData(bcc_cancer)
bcc_cancer <- RunPCA(bcc_cancer)
bcc_cancer <- FindNeighbors(bcc_cancer, dims = 1:20)
bcc_cancer <- RunUMAP(bcc_cancer, dims = 1:20)

# Find optimal pK
sweep.res.bcc <- paramSweep(bcc_cancer, PCs = 1:20, sct = FALSE)
sweep.stats.bcc <- summarizeSweep(sweep.res.bcc, GT = FALSE)
bcmvn.bcc <- find.pK(sweep.stats.bcc)
optimal.pk.bcc <- as.numeric(as.character(bcmvn.bcc$pK[which.max(bcmvn.bcc$BCmetric)]))

# Run DoubletFinder
nExp_poi.bcc <- round(0.028 * ncol(bcc_cancer))
bcc_cancer <- doubletFinder(bcc_cancer, PCs = 1:20, pN = 0.25, pK = optimal.pk.bcc, 
                            nExp = nExp_poi.bcc, sct = FALSE)

# Filter to singlets
df_col_bcc <- grep("DF.classifications", colnames(bcc_cancer@meta.data), value = TRUE)
bcc_singlets <- subset(bcc_cancer, cells = rownames(bcc_cancer@meta.data)[bcc_cancer@meta.data[[df_col_bcc]] == "Singlet"])

# --- LIVER DoubletFinder ---
liver <- JoinLayers(liver)
liver <- NormalizeData(liver)
liver <- FindVariableFeatures(liver, selection.method = "vst", nfeatures = 2000)
liver <- ScaleData(liver)
liver <- RunPCA(liver)
liver <- FindNeighbors(liver, dims = 1:20)
liver <- RunUMAP(liver, dims = 1:20)
# Find optimal pK
sweep.res.liver <- paramSweep(liver, PCs = 1:20, sct = FALSE)
sweep.stats.liver <- summarizeSweep(sweep.res.liver, GT = FALSE)
bcmvn.liver <- find.pK(sweep.stats.liver)
optimal.pk.liver <- as.numeric(as.character(bcmvn.liver$pK[which.max(bcmvn.liver$BCmetric)]))

# Run DoubletFinder
nExp_poi.liver <- round(0.08 * ncol(liver))
liver <- doubletFinder(liver, PCs = 1:20, pN = 0.25, pK = optimal.pk.liver, 
                       nExp = nExp_poi.liver, sct = FALSE)

# Filter to singlets
df_col_liver <- grep("DF.classifications", colnames(liver@meta.data), value = TRUE)
liver_singlets <- subset(liver, cells = rownames(liver@meta.data)[liver@meta.data[[df_col_liver]] == "Singlet"])

# --- LIVER: Filter to Malignant cells ---
samples1 <- read.table(gzfile(paste0(raw_path, "liver/GSE125449_Set1_samples.txt.gz")), header = TRUE, sep = "\t")
samples2 <- read.table(gzfile(paste0(raw_path, "liver/GSE125449_Set2_samples.txt.gz")), header = TRUE, sep = "\t")
all_samples <- rbind(samples1, samples2)

liver_barcodes <- gsub("^(HCC|iCCA)_", "", colnames(liver_singlets))
barcode_to_type <- setNames(all_samples$Type, all_samples$Cell.Barcode)
cell_type_vector <- barcode_to_type[liver_barcodes]
names(cell_type_vector) <- colnames(liver_singlets)
liver_singlets$cell_type <- cell_type_vector

liver_malignant <- subset(liver_singlets, cell_type == "Malignant cell")

# ============================================
# 5. FULL PROCESSING PIPELINE
# ============================================

# --- MELANOMA ---
melanoma_cancer <- NormalizeData(melanoma_cancer)
melanoma_cancer <- FindVariableFeatures(melanoma_cancer, selection.method = "vst", nfeatures = 2000)
melanoma_cancer <- ScaleData(melanoma_cancer)
melanoma_cancer <- RunPCA(melanoma_cancer)
melanoma_cancer <- FindNeighbors(melanoma_cancer, dims = 1:20)
melanoma_cancer <- FindClusters(melanoma_cancer, resolution = 0.5)
melanoma_cancer <- RunUMAP(melanoma_cancer, dims = 1:20)

# --- BCC (already processed, just add clustering) ---
bcc_singlets <- FindClusters(bcc_singlets, resolution = 0.5)

# --- LIVER ---
liver_malignant <- NormalizeData(liver_malignant)
liver_malignant <- FindVariableFeatures(liver_malignant, selection.method = "vst", nfeatures = 2000)
liver_malignant <- ScaleData(liver_malignant)
liver_malignant <- RunPCA(liver_malignant)
liver_malignant <- FindNeighbors(liver_malignant, dims = 1:20)
liver_malignant <- FindClusters(liver_malignant, resolution = 0.5)
liver_malignant <- RunUMAP(liver_malignant, dims = 1:20)

# ============================================
# 6. DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================

# Melanoma: post vs pre treatment
Idents(melanoma_cancer) <- "treatment.group"
melanoma_de <- FindMarkers(melanoma_cancer, ident.1 = "post.treatment", ident.2 = "treatment.naive",
                           min.pct = 0.1, logfc.threshold = 0.25)

# BCC: post vs pre treatment
Idents(bcc_singlets) <- "treatment"
bcc_de <- FindMarkers(bcc_singlets, ident.1 = "post", ident.2 = "pre",
                      min.pct = 0.1, logfc.threshold = 0.25)

# Liver: HCC vs iCCA
Idents(liver_malignant) <- "Cancer_type"
liver_de <- FindMarkers(liver_malignant, ident.1 = "HCC", ident.2 = "iCCA",
                        min.pct = 0.1, logfc.threshold = 0.25)

# Remove mitochondrial genes
bcc_de_filtered <- bcc_de[!grepl("^MT-", rownames(bcc_de)), ]
liver_de_filtered <- liver_de[!grepl("^MT-|^MTRNR", rownames(liver_de)), ]

# ============================================
# 7. VISUALIZATIONS
# ============================================

# UMAP plots
p1_clusters <- DimPlot(melanoma_cancer, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("Melanoma - Clusters")
p1_treatment <- DimPlot(melanoma_cancer, reduction = "umap", group.by = "treatment.group") + ggtitle("Melanoma - Treatment")
melanoma_umap <- p1_clusters + p1_treatment
ggsave(paste0(fig_path, "melanoma_umap.png"), melanoma_umap, width = 12, height = 5)

p2_clusters <- DimPlot(bcc_singlets, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("BCC - Clusters")
p2_treatment <- DimPlot(bcc_singlets, reduction = "umap", group.by = "treatment") + ggtitle("BCC - Treatment")
bcc_umap <- p2_clusters + p2_treatment
ggsave(paste0(fig_path, "bcc_umap.png"), bcc_umap, width = 12, height = 5)

p3_clusters <- DimPlot(liver_malignant, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("Liver - Clusters")
p3_cancertype <- DimPlot(liver_malignant, reduction = "umap", group.by = "Cancer_type") + ggtitle("Liver - Cancer Type")
liver_umap <- p3_clusters + p3_cancertype
ggsave(paste0(fig_path, "liver_umap.png"), liver_umap, width = 12, height = 5)

# Volcano plots
create_volcano <- function(de_data, title) {
  de_data$gene <- rownames(de_data)
  de_data$significance <- "NS"
  de_data$significance[de_data$avg_log2FC > 0.5 & de_data$p_val_adj < 0.05] <- "Up"
  de_data$significance[de_data$avg_log2FC < -0.5 & de_data$p_val_adj < 0.05] <- "Down"
  
  ggplot(de_data, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("Down" = "blue", "NS" = "grey", "Up" = "red")) +
    theme_minimal() +
    ggtitle(title) +
    xlab("Log2 Fold Change") +
    ylab("-Log10 Adjusted P-value")
}
v1 <- create_volcano(melanoma_de, "Melanoma: Post vs Pre Treatment")
v2 <- create_volcano(bcc_de_filtered, "BCC: Post vs Pre Treatment")
v3 <- create_volcano(liver_de_filtered, "Liver: HCC vs iCCA")

ggsave(paste0(fig_path, "volcano_melanoma.png"), v1, width = 8, height = 6)
ggsave(paste0(fig_path, "volcano_bcc.png"), v2, width = 8, height = 6)
ggsave(paste0(fig_path, "volcano_liver.png"), v3, width = 8, height = 6)

# ============================================
# 8. SAVE RESULTS
# ============================================

# Save processed objects
saveRDS(melanoma_cancer, paste0(processed_path, "melanoma_final_processed.RDS"))
saveRDS(bcc_singlets, paste0(processed_path, "bcc_final_processed.RDS"))
saveRDS(liver_malignant, paste0(processed_path, "liver_final_processed.RDS"))

# Save DE tables
write.csv(melanoma_de, paste0(table_path, "melanoma_DE_post_vs_pre.csv"))
write.csv(bcc_de_filtered, paste0(table_path, "bcc_DE_post_vs_pre_noMT.csv"))
write.csv(liver_de_filtered, paste0(table_path, "liver_DE_HCC_vs_iCCA_noMT.csv"))
# ============================================
# 9. FINAL SUMMARY
# ============================================

print("========================================")
print("   ICB RAW REPRODUCTION - COMPLETE!    ")
print("========================================")
print(paste("Melanoma cells:", ncol(melanoma_cancer), "(Paper: 1,945)"))
print(paste("BCC cells:", ncol(bcc_singlets), "(Paper: 3,500)"))
print(paste("Liver cells:", ncol(liver_malignant), "(Paper: 1,992)"))
print("========================================")

