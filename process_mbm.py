#!/usr/bin/env python

"""
Single-cell RNA-seq analysis pipeline for Alvarez-Breckenridge MBM data set
Python/Scanpy conversion of the original R/Seurat pipeline
"""
# Import necessary libraries
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def main():
    # --- Set up Scanpy ---
    sc.settings.verbosity = 3  # Show detailed logs
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=80, facecolor='white')

    print("=" * 80)
    print("MBM Dataset Analysis Pipeline")
    print("=" * 80)

    # ============================================================================
    # --- 1. Reading the scRNA-seq data ---
    # ============================================================================
    print("--- Reading the scRNA-seq data ---")

    # Set the path to "./" to read from the current directory
    path_to_data = "./"

    # Define file paths
    meta_file1 = f"{path_to_data}BM_c1000_base_metadata.tsv"
    mtx_file = f"{path_to_data}BM_c1000_raw.mtx"
    genes_file = f"{path_to_data}BM_c1000_preQC_genes_raw.tsv"
    cells_file = f"{path_to_data}BM_c1000_preQC_cellnames_raw.tsv"
    meta_file2 = f"{path_to_data}postQC_c1000_spatial.tsv"

    # Set the exact name of your Excel file
    meta_file3_excel = f"{path_to_data}cir-21-0870_supplementary_data_3_suppsm3.xlsx"

    print(f"Loading Excel file: {meta_file3_excel}...")

    # Use pd.read_excel() to read the file.
    clinical_data_raw = pd.read_excel(meta_file3_excel, header=None)

    # Transpose the data
    clinical_data = clinical_data_raw.T
    # Set the first row as the new header
    clinical_data.columns = clinical_data.iloc[0]
    # Get the patient IDs from the index
    clinical_data['donor_id'] = clinical_data.index
    # Remove the first row (which is now the header)
    clinical_data = clinical_data.iloc[1:]

    # Load first metadata file
    print(f"Loading {meta_file1}...")
    base_metadata = pd.read_csv(
        meta_file1,
        sep='\t'  # 1. Read normally (header is first row)
    )
    base_metadata = base_metadata.iloc[1:]  # 2. Remove first data row (R: [-c(1),])
    base_metadata = base_metadata.set_index('NAME')  # 3. Set index from 'NAME' column

    # Load the raw count matrix
    print(f"Loading {mtx_file}...")
    counts = sc.read_mtx(mtx_file).T

    # Load cell barcodes (obs_names)
    print(f"Loading {cells_file}...")
    barcodes = pd.read_csv(cells_file, sep='\t', header=None)[0].values

    # Load gene names (var_names)
    print(f"Loading {genes_file}...")
    genes = pd.read_csv(genes_file, sep='\t', header=None, usecols=[0])[0].values

    # Create the AnnData object
    print("Creating AnnData object...")
    adata = sc.AnnData(counts)
    adata.obs_names = barcodes
    adata.var_names = genes

    # Add the first metadata
    # We use .join, as .loc might fail if indices are not perfectly aligned
    adata.obs = adata.obs.join(base_metadata)

    # Load and add the second metadata file (clustering)
    print(f"Loading {meta_file2}...")
    clustering_data = pd.read_csv(
        meta_file2,
        sep='\t'  # 1. Read normally
    )
    clustering_data = clustering_data.iloc[1:]  # 2. Remove first data row
    clustering_data = clustering_data.set_index('NAME')  # 3. Set index

    adata.obs = adata.obs.join(clustering_data)

    # Merge clinical data into the AnnData object
    adata.obs = adata.obs.reset_index().rename(columns={'index': 'cell_id'})
    # --- START OF NEW FIX ---
    # Force both 'donor_id' columns to be strings (object)
    print("Fixing donor_id data types for merge...")
    if 'donor_id' in adata.obs.columns:
        adata.obs['donor_id'] = adata.obs['donor_id'].astype(str)
        clinical_data['donor_id'] = clinical_data['donor_id'].astype(str)
    else:
        print("Error: 'donor_id' column not found in base metadata.")
    # --- END OF NEW FIX ---
    adata.obs = adata.obs.merge(clinical_data, on='donor_id', how='left')
    adata.obs.index = adata.obs['cell_id']

    print("Raw data loaded and merged:")
    print(adata)

    # --- 2. Refining the seurat object ---
    print("--- 2. Refining the seurat object ---")

    # Clean the 'cell type' column
    adata.obs['cell type'] = adata.obs['cell type'].fillna('Unknown')
    adata.obs['cell type'] = adata.obs['cell type'].replace('', 'Unknown')

    # Filter out "Unknown" cells
    adata = adata[adata.obs['cell type'] != 'Unknown', :].copy()

    # Rename column
    adata.obs = adata.obs.rename(columns={'cell type': 'cell_types'})

    # Standardize the 'outcome' column
    if 'outcome' in adata.obs.columns:
        adata.obs['outcome'] = adata.obs['outcome'].replace('No ICI', 'UT')
    else:
        # Check for the column name from your output, which has a space
        outcome_col = 'Responder (R), Partial respodner (PR),  Non-responnder (NR), No ICI '
        if outcome_col in adata.obs.columns:
            print(f"Found outcome column: '{outcome_col}'")
            adata.obs[outcome_col] = adata.obs[outcome_col].replace('No ICI', 'UT')
            # Rename it to something simple
            adata.obs = adata.obs.rename(columns={outcome_col: 'outcome'})
        else:
            print("Warning: 'outcome' column not found. Skipping outcome standardization.")

    # Create the 'pre_post' column
    if 'Pre/post ICI' in adata.obs.columns:
        adata.obs['pre_post'] = np.where(adata.obs['Pre/post ICI'] == 'Post', 'Post', 'Pre')
    else:
        print("Warning: 'Pre/post ICI' column not found. Cannot create 'pre_post' column.")

    # Create combined metadata columns
    adata.obs['donor_id_pre_post'] = adata.obs['donor_id'] + '_' + adata.obs['pre_post']
    adata.obs['donor_id_outcome'] = adata.obs['donor_id'] + '_' + adata.obs['outcome'].astype(str)
    adata.obs['donor_id_cell_types'] = adata.obs['donor_id'] + '_' + adata.obs['cell_types']
    adata.obs['donor_id_cell_types_pre_post'] = adata.obs['donor_id_cell_types'] + '_' + adata.obs['pre_post']
    adata.obs['sample_id_pre_post_outcome'] = adata.obs['donor_id_pre_post'] + '_' + adata.obs['outcome'].astype(str)

    # --- R: removing samples with less than 20 cells ---
    print("Applying > 20 cells per group filter...")

    # 1. Get counts for each group
    group_counts = adata.obs.groupby('donor_id_cell_types_pre_post')['cell_id'].count()
    # 2. Find groups with > 20 cells
    groups_to_keep = group_counts[group_counts > 20].index
    # 3. Create the 'enough_cells' column
    adata.obs['enough_cells'] = np.where(
        adata.obs['donor_id_cell_types_pre_post'].isin(groups_to_keep),
        'enough',
        'not_enough'
    )

    # Filter the AnnData object. .copy() is important!
    adata_subset = adata[adata.obs['enough_cells'] == 'enough', :].copy()

    # Add standardized study information
    adata_subset.obs['Study_name'] = 'Alvarez_Breckenridge'
    adata_subset.obs['Cancer_type'] = 'Melanoma_derived_brain_metastases'
    adata_subset.obs['Primary_or_met'] = 'Metastatic'

    print(f"Original cell count: {adata.n_obs}")
    print(f"Filtered cell count (> 20 per group): {adata_subset.n_obs}")

    # --- 3. scRNA-seq analysis (clustering + normalization) ---
    print("--- 3. scRNA-seq analysis (clustering + normalization) ---")

    # --- R: seurat.object <- NormalizeData(seurat.object) ---
    print("Normalizing and Log-transforming data...")
    sc.pp.normalize_total(adata_subset, target_sum=1e4)
    sc.pp.log1p(adata_subset)

    # --- R: seurat.object <- FindVariableFeatures(seurat.object) ---
    print("Finding highly variable genes...")
    sc.pp.highly_variable_genes(adata_subset, min_mean=0.0125, max_mean=3, min_disp=0.5)

    # --- R: seurat.object <- ScaleData(seurat.object, ...) ---
    print("Scaling data...")
    sc.pp.scale(adata_subset, max_value=10)

    # --- R: seurat.object <- RunPCA(seurat.object) ---
    print("Running PCA...")
    sc.tl.pca(adata_subset, svd_solver='arpack')

    # --- R: ElbowPlot(seurat.object) ---
    # We'll save this plot to the current folder
    print("Running PCA and saving variance ratio plot (linear scale)...")
    # 💡 MODIFICATION: log=False for a linear scale, which shows the 'elbow'
    elbow_plot = sc.pl.pca_variance_ratio(adata_subset, log=False, show=False)
    plt.savefig('pca_variance_ratio_linear.png')
    plt.close()  # Close the plot to save memory
    print("Saved pca_variance_ratio_linear.png")

    # --- R: seurat.object <- FindNeighbors(seurat.object, dims = 1:20) ---
    print("Computing neighbors graph...")
    sc.pp.neighbors(adata_subset, n_neighbors=10, n_pcs=20)

    # --- R: seurat.object <- FindClusters(seurat.object) ---
    print("Finding clusters (Leiden)...")
    sc.tl.leiden(adata_subset)

    # --- R: seurat.object <- RunUMAP(seurat.object, dims = 1:20) ---
    print("Running UMAP...")
    sc.tl.umap(adata_subset)

    print("--- Analysis Complete ---")
    print("Final processed object:")
    print(adata_subset)

    # Convert all columns that might contain mixed types/IDs to explicit strings
    for col in adata_subset.obs.columns:
        if adata_subset.obs[col].dtype == 'object':
            # This is slow, but fixes the H5AD saving issue on complex metadata
            adata_subset.obs[col] = adata_subset.obs[col].astype(str)

    print("Saving UMAP plots: ")

    # Plot 1: Leiden Clusters (for identifying groups)
    sc.pl.umap(
        adata_subset,
        color='leiden',
        title='Leiden Clustering',
        save='_leiden.png',
        show=False,
        # Adjust aesthetics for large dataset visualization
        s=15,  # Marker size (smaller dots)
        alpha=0.6,  # Transparency (helps see density)
        legend_fontsize=8,
        frameon=False  # Removes the border around the plot
    )
    plt.close()

    # Plot 2: Cell Types (for biological interpretation)
    print("Saving UMAP plot for Cell Types with fixed legend location...")
    sc.pl.umap(
        adata_subset,
        color='cell_types',
        title='Cell Type Annotation',
        save='_cell_types_revised.png',  # Saving to a new file name
        show=False,
        s=30,  # Slightly larger marker size (was 15)
        alpha=0.6,
        legend_fontsize=8,
        frameon=False
    )
    plt.close()

    # Plot 3: Pre/Post Status (for experimental factor visualization)
    sc.pl.umap(
        adata_subset,
        color='pre_post',
        title='Pre/Post ICI Status',
        save='_pre_post.png',
        show=False,
        s=15,
        alpha=0.6,
        frameon=False
    )
    plt.close()

    print("Saved UMAP plots: umap_leiden.png, umap_cell_types.png, umap_pre_post.png")

    # --- R: saveRDS(seurat_MBM_Christopher_m_epi_subset, ...) ---
    # Save the final, processed object to an .h5ad file
    output_file = "MBM_Christopher.h5ad"
    print(f"Saving final object to {output_file}...")
    adata_subset.write_h5ad(output_file)

    print(f"Successfully converted and processed. Final object saved to: {output_file}")

if __name__ == "__main__":
    main()
