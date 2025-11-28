# ============================================================================
# tme_mal.R
# - Malignant cell state scoring (manual + bulk-derived metaprograms)
# - TME composition
# ============================================================================

# Load packages
library(Seurat)
library(dplyr)
library(ggplot2)

# Load scRNA objects
mal <- readRDS("Malignant_cells.rds")
tme <- readRDS("Entire_TME.rds")

# Inspect
mal
colnames(mal@meta.data)
head(mal@meta.data)

tme
colnames(tme@meta.data)
head(tme@meta.data)

# Malignant UMAPs by existing metadata
DimPlot(mal, reduction = "umap", group.by = "Malignant_clusters",
        label = TRUE, repel = TRUE) +
  ggtitle("Malignant clusters")

DimPlot(mal, reduction = "umap", group.by = "Response") +
  ggtitle("Response")

DimPlot(mal, reduction = "umap", group.by = "Timepoint") +
  ggtitle("Timepoint")

# TME UMAPs by cluster / sample / timepoint

DimPlot(tme, reduction = "umap", group.by = "seurat_clusters",
        label = TRUE, repel = TRUE) +
  ggtitle("TME clusters")

DimPlot(tme, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Sample ID")

DimPlot(tme, reduction = "umap", group.by = "BT/OT") +
  ggtitle("Timepoint / Treatment (BT vs OT)")

# Malignant cells - MANUAL metaprograms (MEL, MES, AP, NC, Stress)

manual_metaprograms <- list(
  MEL       = c("MITF", "PMEL", "TYR"),
  MES       = c("AXL", "NGFR", "WNT5A"),
  AP        = c("CDK1", "TOP2A", "MKI67"),
  NC        = c("SOX10", "EDNRB"),
  StressH   = c("HSPA1A", "HSP90AA1"),
  StressP53 = c("TP53", "CDKN1A")
)

for (prog in names(manual_metaprograms)) {
  mal <- AddModuleScore(
    mal,
    features = list(manual_metaprograms[[prog]]),
    name = paste0(prog, "_score")
  )
}

# Identify all manual score columns (e.g., "MEL_score1", "MES_score1", ...)
manual_score_cols <- grep("_score", colnames(mal@meta.data), value = TRUE)

# Dominant manual state per cell
mal$dominant_manual_state <- apply(mal@meta.data[, manual_score_cols, drop = FALSE], 1, function(x) {
  names(x)[which.max(x)]
})

DimPlot(mal, reduction = "umap", group.by = "dominant_manual_state",
        label = TRUE, repel = TRUE) +
  ggtitle("Malignant cell metaprograms (manual signatures)")

# Fractions per patient (manual states)
mal_manual_meta <- mal@meta.data %>%
  group_by(patient_ID, dominant_manual_state, Response) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(patient_ID) %>%
  mutate(fraction = n / sum(n))

ggplot(mal_manual_meta,
       aes(x = patient_ID, y = fraction, fill = dominant_manual_state)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Response, scales = "free_x") +
  theme_minimal() +
  ylab("Fraction of malignant cells") +
  xlab("Patient ID") +
  ggtitle("Manual metaprogram fractions per patient") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Malignant cells - BULK-DERIVED metaprograms from NMF

# This file is created by bulk_metaprograms.R
if (file.exists("Bulk_metaprogram_genes_top50.rds")) {
  bulk_metaprograms <- readRDS("Bulk_metaprogram_genes_top50.rds")
  
  # Score each bulk metaprogram
  mal <- AddModuleScore(
    mal,
    features = bulk_metaprograms,
    name = "BulkProg_"
  )
  
  bulk_score_cols <- grep("^BulkProg_", colnames(mal@meta.data), value = TRUE)
  
  # Dominant bulk program per cell
  mal$bulk_dominant_prog <- apply(mal@meta.data[, bulk_score_cols, drop = FALSE], 1, function(x) {
    names(x)[which.max(x)]
  })
  mal$bulk_dominant_prog <- gsub("BulkProg_", "Prog", mal$bulk_dominant_prog)
  
  DimPlot(mal, reduction = "umap", group.by = "bulk_dominant_prog",
          label = TRUE, repel = TRUE) +
    ggtitle("Malignant cells – dominant bulk-derived metaprogram")
  
  # Fractions per patient (bulk programs)
  mal_bulk_meta <- mal@meta.data %>%
    group_by(patient_ID, bulk_dominant_prog, Response) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(patient_ID) %>%
    mutate(fraction = n / sum(n))
  
  ggplot(mal_bulk_meta,
         aes(x = patient_ID, y = fraction, fill = bulk_dominant_prog)) +
    geom_bar(stat = "identity") +
    facet_wrap(~Response, scales = "free_x") +
    theme_minimal() +
    ylab("Fraction of malignant cells") +
    xlab("Patient ID") +
    ggtitle("Bulk-derived metaprogram fractions per patient") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# TME annotation (cluster → cell type mapping)

cluster2celltype <- c(
  "0" = "T cells",
  "1" = "Macrophages",
  "2" = "B cells",
  "3" = "NK cells",
  "4" = "Fibroblasts",
  "5" = "Endothelial cells",
  "6" = "Dendritic cells"
)

tme$celltype <- unname(cluster2celltype[as.character(tme$seurat_clusters)])

DimPlot(tme, reduction = "umap", group.by = "celltype",
        label = TRUE, repel = TRUE) +
  ggtitle("TME cell types")

# Optional marker validation
VlnPlot(tme,
        features = c("CD3D","MS4A1","CD68","GNLY","COL1A1","PECAM1"),
        group.by = "celltype", pt.size = 0)

# TME composition per sample / BT vs OT

tme_meta <- tme@meta.data %>%
  group_by(orig.ident, celltype, `BT/OT`) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(fraction = n / sum(n))

ggplot(tme_meta,
       aes(x = orig.ident, y = fraction, fill = celltype)) +
  geom_bar(stat = "identity") +
  facet_wrap(~`BT/OT`, scales = "free_x") +
  theme_minimal() +
  ylab("Fraction of TME cells") +
  xlab("Sample ID") +
  ggtitle("TME composition per sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save prepared objects

saveRDS(mal, file = "mal_prepared_all.rds")   # manual + bulk scores
saveRDS(tme, file = "tme_prepared_step3.rds")
