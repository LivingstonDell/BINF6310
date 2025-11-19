# ==========================================================
# Integrated Workflow: Step 2 + Step 3
# Malignant and TME cells
# ==========================================================

# Load packages
library(Seurat)
library(dplyr)
library(ggplot2)

# Load objects
mal <- readRDS("Malignant_cells.rds")
tme <- readRDS("Entire_TME.rds")

# Inspect Malignant Cells
mal
colnames(mal@meta.data)
head(mal@meta.data)
mal@reductions

# Visualize UMAP by existing metadata
DimPlot(mal, reduction = "umap", group.by = "Malignant_clusters", label = TRUE, repel = TRUE) +
  ggtitle("Malignant clusters")

DimPlot(mal, reduction = "umap", group.by = "Response") +
  ggtitle("Response")

DimPlot(mal, reduction = "umap", group.by = "Timepoint") +
  ggtitle("Timepoint")

# Inspect TME
tme
colnames(tme@meta.data)
head(tme@meta.data)
tme@reductions

# Visualize UMAP by clusters and sample/timepoint
DimPlot(tme, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  ggtitle("TME clusters")

DimPlot(tme, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Sample ID")

DimPlot(tme, reduction = "umap", group.by = "BT/OT") +
  ggtitle("Timepoint / Treatment")


# Malignant cells - assign metaprograms

# Define marker genes
metaprograms <- list(
  MEL = c("MITF", "PMEL", "TYR"),
  MES = c("AXL", "NGFR", "WNT5A"),
  AP  = c("CDK1", "TOP2A", "MKI67"),
  NC  = c("SOX10", "EDNRB"),
  StressH = c("HSPA1A", "HSP90AA1"),
  StressP53 = c("TP53", "CDKN1A")
)

# Score each metaprogram
for (prog in names(metaprograms)) {
  mal <- AddModuleScore(
    mal,
    features = list(metaprograms[[prog]]),
    name = paste0(prog, "_score")
  )
}

# Assign dominant metaprogram per cell
score_cols <- grep("_score", colnames(mal@meta.data), value = TRUE)
mal$dominant_state <- apply(mal@meta.data[, score_cols], 1, function(x) {
  names(x)[which.max(x)]
})

# UMAP plot by dominant_state
DimPlot(mal, reduction = "umap", group.by = "dominant_state", label = TRUE, repel = TRUE) +
  ggtitle("Malignant cell metaprograms")

# Compute fractions per patient
mal_meta <- mal@meta.data %>%
  group_by(patient_ID, dominant_state, Response) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(patient_ID) %>%
  mutate(fraction = n / sum(n))

# Stacked barplot
ggplot(mal_meta, aes(x = patient_ID, y = fraction, fill = dominant_state)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Response, scales = "free_x") +
  theme_minimal() +
  ylab("Fraction of malignant cells") +
  xlab("Patient ID") +
  ggtitle("Metaprogram fractions per patient")

# Step 3b: TME annotation
# Map clusters to canonical cell types (adjust based on your clusters)
cluster2celltype <- c(
  "0" = "T cells",
  "1" = "Macrophages",
  "2" = "B cells",
  "3" = "NK cells",
  "4" = "Fibroblasts",
  "5" = "Endothelial cells",
  "6" = "Dendritic cells"
)
tme$celltype <- cluster2celltype[as.character(tme$seurat_clusters)]

# UMAP colored by cell type
DimPlot(tme, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) +
  ggtitle("TME cell types")

# Optional: violin plot to validate markers
VlnPlot(tme, features = c("CD3D","MS4A1","CD68","GNLY","COL1A1","PECAM1"),
        group.by = "celltype", pt.size = 0)

# Compute fractions per sample
tme_meta <- tme@meta.data %>%
  group_by(orig.ident, celltype, `BT/OT`) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(fraction = n / sum(n))

# Stacked barplot of TME composition
ggplot(tme_meta, aes(x = orig.ident, y = fraction, fill = celltype)) +
  geom_bar(stat = "identity") +
  facet_wrap(~`BT/OT`, scales = "free_x") +
  theme_minimal() +
  ylab("Fraction of TME cells") +
  xlab("Sample ID") +
  ggtitle("TME composition per sample")

# Save prepared objects
saveRDS(mal, file = "mal_prepared_step3.rds")
saveRDS(tme, file = "tme_prepared_step3.rds")
