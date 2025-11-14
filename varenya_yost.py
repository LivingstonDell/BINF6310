#!/usr/bin/env python3
"""
Single-cell RNA-seq analysis pipeline for BCC/SCC Yost dataset
Python/Scanpy conversion of the original R/Seurat pipeline
"""

import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
from scipy.sparse import csr_matrix
import warnings
warnings.filterwarnings('ignore')

# Set scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

print("=" * 80)
print("BCC/SCC Yost Dataset Analysis Pipeline")
print("=" * 80)

# ============================================================================
# 1. READING THE scRNA-seq DATA
# ============================================================================
print("\n[1/6] Reading BCC data...")

# Read BCC metadata
bcc_metadata = pd.read_csv(
    "GSE123813_bcc_all_metadata.txt.gz",
    sep="\t",
    compression='gzip'
)
bcc_metadata.set_index('cell.id', inplace=True)

# Read BCC counts
bcc_counts = pd.read_csv(
    "GSE123813_bcc_scRNA_counts.txt.gz",
    sep="\t",
    compression='gzip',
    index_col=0
)

# Create AnnData object for BCC
adata_bcc = ad.AnnData(
    X=csr_matrix(bcc_counts.T),  # Transpose: cells x genes
    obs=bcc_metadata,
    var=pd.DataFrame(index=bcc_counts.index)
)
adata_bcc.obs['dataset'] = 'BCC'

print(f"  BCC: {adata_bcc.n_obs} cells, {adata_bcc.n_vars} genes")

# Read SCC data
print("\n[2/6] Reading SCC data...")

scc_metadata = pd.read_csv(
    "GSE123813_scc_metadata.txt.gz",
    sep="\t",
    compression='gzip'
)
scc_metadata.set_index('cell.id', inplace=True)

scc_counts = pd.read_csv(
    "GSE123813_scc_scRNA_counts.txt.gz",
    sep="\t",
    compression='gzip',
    index_col=0
)

# Create AnnData object for SCC
adata_scc = ad.AnnData(
    X=csr_matrix(scc_counts.T),
    obs=scc_metadata,
    var=pd.DataFrame(index=scc_counts.index)
)
adata_scc.obs['dataset'] = 'SCC'

print(f"  SCC: {adata_scc.n_obs} cells, {adata_scc.n_vars} genes")

# Merge BCC and SCC datasets
print("\n[3/6] Merging BCC and SCC datasets...")
adata_combined = ad.concat([adata_bcc, adata_scc], join='outer', fill_value=0)
print(f"  Combined: {adata_combined.n_obs} cells, {adata_combined.n_vars} genes")

# ============================================================================
# 2. CALCULATING CELL TYPE COUNTS AND IMMUNE INFILTRATION METRICS
# ============================================================================
print("\n[4/6] Calculating cell type counts and immune metrics...")

# Define CD8 T cell types
cd8_types = [
    "CD8_act", "CD8_act_T_cells", "CD8_eff", "CD8_ex", 
    "CD8_ex_act", "CD8_ex_T_cells", "CD8_mem", "CD8_mem_T_cells", "CD8_naive"
]

# Define all T cell types
t_cell_types = cd8_types + [
    "CD4_T_cells", "Tcell_prolif", "Tfh", "Th17", "Treg", "Tregs"
]

# Create unified CD8 label
adata_combined.obs['cluster_cd8'] = adata_combined.obs['cluster'].apply(
    lambda x: 'CD8t' if x in cd8_types else x
)

# Create unified T cell label
adata_combined.obs['cluster_tcell'] = adata_combined.obs['cluster'].apply(
    lambda x: 'T_cells' if x in t_cell_types else x
)

# Calculate CD8 counts per patient
cd8_counts = adata_combined.obs[adata_combined.obs['cluster_cd8'] == 'CD8t'].groupby('patient').size()
cd8_counts = cd8_counts.to_frame('total_T_Cell').reset_index()
cd8_counts.columns = ['sample_id', 'total_T_Cell']

# Calculate T cell counts per patient
tcell_counts = adata_combined.obs[adata_combined.obs['cluster_tcell'] == 'T_cells'].groupby('patient').size()
tcell_counts = tcell_counts.to_frame('total_T_Cell_only').reset_index()
tcell_counts.columns = ['sample_id', 'total_T_Cell_only']

# Merge CD8 and T cell counts
immune_counts = pd.merge(cd8_counts, tcell_counts, on='sample_id', how='outer')

# Calculate total cells per patient
total_counts = adata_combined.obs.groupby('patient').size()
total_counts = total_counts.to_frame('total_cell_per_patient').reset_index()
total_counts.columns = ['sample_id', 'total_cell_per_patient']

# Merge all counts
cell_counts = pd.merge(total_counts, immune_counts, on='sample_id', how='left')
cell_counts.fillna(0, inplace=True)

# ============================================================================
# 3. SUBSETTING TO TUMOR CELLS
# ============================================================================
print("\n[5/6] Subsetting to tumor cells and filtering...")

# Identify epithelial/tumor cells
adata_combined.obs['epi'] = adata_combined.obs['cluster'].isin(['Tumor_1', 'Tumor_2'])
adata_tumor = adata_combined[adata_combined.obs['epi']].copy()

print(f"  Tumor cells: {adata_tumor.n_obs}")

# Filter: keep only patients with >20 tumor cells
patient_counts = adata_tumor.obs['patient'].value_counts()
patients_enough = patient_counts[patient_counts > 20].index
adata_tumor = adata_tumor[adata_tumor.obs['patient'].isin(patients_enough)].copy()

print(f"  After filtering (>20 cells/patient): {adata_tumor.n_obs} cells")

# Relabel treatment status
adata_tumor.obs['pre_post'] = adata_tumor.obs['treatment'].map({'pre': 'Pre', 'post': 'Post'})

# Add study name
adata_tumor.obs['Study_name'] = 'Yost'
adata_tumor.obs['Primary_or_met'] = 'Metastatic'

# Rename patient column to sample_id
adata_tumor.obs['sample_id'] = adata_tumor.obs['patient']

# Create sample_id with timepoint
adata_tumor.obs['sample_id_pre_post'] = (
    adata_tumor.obs['sample_id'] + '_' + adata_tumor.obs['pre_post']
)

# Create outcome column (without clinical response data, only timepoint info)
adata_tumor.obs['outcome'] = adata_tumor.obs['pre_post'].apply(
    lambda x: 'UT' if x == 'Pre' else 'Post'  # UT = Untreated
)

adata_tumor.obs['sample_id_outcome'] = (
    adata_tumor.obs['sample_id'] + '_' + adata_tumor.obs['outcome']
)

# Filter: keep only sample_id_pre_post with >20 cells
sample_counts = adata_tumor.obs['sample_id_pre_post'].value_counts()
samples_enough = sample_counts[sample_counts > 20].index
adata_tumor = adata_tumor[adata_tumor.obs['sample_id_pre_post'].isin(samples_enough)].copy()

print(f"  After filtering (>20 cells/sample_timepoint): {adata_tumor.n_obs} cells")

# Rename cell.id to cell_id for consistency
adata_tumor.obs['cell_id'] = adata_tumor.obs.index

# ============================================================================
# 4. ADDING IMMUNE CELL METRICS
# ============================================================================
print("\n  Adding immune infiltration metrics...")

# Filter cell_counts to only include samples in final dataset
samples_in_data = adata_tumor.obs['sample_id'].unique()
cell_counts_filtered = cell_counts[cell_counts['sample_id'].isin(samples_in_data)].copy()

# Merge immune metrics into adata
adata_tumor.obs = adata_tumor.obs.merge(
    cell_counts_filtered,
    on='sample_id',
    how='left'
)

# Calculate normalized CD8 abundance (CD8/total cells)
adata_tumor.obs['normalized_CD8_actual_totalcells'] = (
    adata_tumor.obs['total_T_Cell'] / adata_tumor.obs['total_cell_per_patient']
)

# Calculate normalized T cell abundance (T cells/total cells)
adata_tumor.obs['normalized_CD8_totalcells'] = (
    adata_tumor.obs['total_T_Cell_only'] / adata_tumor.obs['total_cell_per_patient']
)

# ============================================================================
# 5. NORMALIZATION AND DIMENSIONALITY REDUCTION
# ============================================================================
print("\n[6/6] Running normalization and dimensionality reduction...")

# Normalize data (log1p transformation, like Seurat's NormalizeData)
sc.pp.normalize_total(adata_tumor, target_sum=1e4)
sc.pp.log1p(adata_tumor)

# Store raw counts
adata_tumor.raw = adata_tumor

# Find highly variable genes
sc.pp.highly_variable_genes(adata_tumor, n_top_genes=2000)

# Scale data
sc.pp.scale(adata_tumor, max_value=10)

# PCA
sc.tl.pca(adata_tumor, svd_solver='arpack')

# Neighbors and clustering
sc.pp.neighbors(adata_tumor, n_neighbors=10, n_pcs=20)
sc.tl.leiden(adata_tumor, resolution=0.8)

# UMAP
sc.tl.umap(adata_tumor, min_dist=0.3)

print(f"\n  PCA computed: {adata_tumor.obsm['X_pca'].shape[1]} components")
print(f"  UMAP computed: {adata_tumor.obsm['X_umap'].shape}")
print(f"  Clusters identified: {adata_tumor.obs['leiden'].nunique()}")

# ============================================================================
# 6. SAVE OUTPUT
# ============================================================================
output_file = "seurat_BCC_SCC_tumor_subset_subset.h5ad"
print(f"\n{'=' * 80}")
print(f"Saving processed data to: {output_file}")
adata_tumor.write(output_file)

# Print summary statistics
print(f"\n{'=' * 80}")
print("ANALYSIS SUMMARY")
print(f"{'=' * 80}")
print(f"Total cells: {adata_tumor.n_obs}")
print(f"Total genes: {adata_tumor.n_vars}")
print(f"Number of patients: {adata_tumor.obs['sample_id'].nunique()}")
print(f"Number of samples (patient_timepoint): {adata_tumor.obs['sample_id_pre_post'].nunique()}")
print(f"\nSample distribution:")
print(adata_tumor.obs.groupby(['sample_id', 'pre_post', 'outcome']).size().to_string())
print(f"\n{'=' * 80}")
print("Pipeline completed successfully!")
print(f"{'=' * 80}\n")

# Optional: Create a simple visualization
print("Generating QC plots...")
sc.pl.umap(adata_tumor, color=['leiden', 'pre_post', 'sample_id'], 
           save='_overview.png', show=False)
print("  Saved: figures/umap_overview.png")