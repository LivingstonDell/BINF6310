import scanpy as sc
import pandas as pd
import os

# Define the file paths
H5AD_FILE = "MBM_Christopher.h5ad"
METADATA_OUTPUT_CSV = "processed_metadata_obs.csv"
NORMALIZED_MATRIX_OUTPUT_TSV = "normalized_expression_matrix.tsv.gz"  # .gz compresses the large matrix file


def export_data():
    """
    Loads the AnnData object and exports its key components (metadata and
    normalized expression matrix) to separate, easy-to-read text files.
    """
    if not os.path.exists(H5AD_FILE):
        print(f"Error: The input file '{H5AD_FILE}' was not found.")
        print("Please ensure the script is run in the same directory as the .h5ad file.")
        return

    print(f"Loading AnnData object from {H5AD_FILE}...")
    try:
        # Load the AnnData object
        adata = sc.read_h5ad(H5AD_FILE)
        print("Successfully loaded AnnData object.")

        # --- 1. Export Metadata (Cell Annotations) ---
        print(f"Exporting metadata to {METADATA_OUTPUT_CSV}...")

        # Pull UMAP coordinates from obsm and add them as columns to the obs dataframe
        # Check if X_umap exists before proceeding
        if 'X_umap' not in adata.obsm:
            print("Warning: X_umap not found in AnnData object. Skipping UMAP export.")
            metadata_df = adata.obs
        else:
            umap_df = pd.DataFrame(adata.obsm['X_umap'],
                                   index=adata.obs.index,
                                   columns=['UMAP_1', 'UMAP_2'])
            # Concatenate the main obs data with the UMAP coordinates
            metadata_df = pd.concat([adata.obs, umap_df], axis=1)

        # Write the combined metadata to a standard CSV file
        metadata_df.to_csv(METADATA_OUTPUT_CSV)
        print(f"Successfully saved metadata (including UMAP and Leiden clusters) to: {METADATA_OUTPUT_CSV}")

        # --- 2. Export Normalized Expression Matrix ---
        print(f"Exporting normalized expression matrix to {NORMALIZED_MATRIX_OUTPUT_TSV}...")

        # FIX: The .h5ad file likely stores X as a dense array, so we remove .todense()

        # We need to handle both possible array types: sparse and dense
        if hasattr(adata.X, 'todense'):
            matrix_data = adata.X.todense()
        else:
            # Assume it's already a dense numpy array
            matrix_data = adata.X

        # Convert the expression matrix to a DataFrame
        # We use .T (transpose) so genes are columns (standard format)
        # We use .tsv.gz for compression, which is crucial for large matrices
        normalized_df = pd.DataFrame(
            matrix_data.T,
            index=adata.var_names,
            columns=adata.obs_names
        ).T  # Final transpose to get Cells as rows, Genes as columns

        normalized_df.to_csv(
            NORMALIZED_MATRIX_OUTPUT_TSV,
            sep='\t',
            compression='gzip'
        )
        print(f"Successfully saved normalized matrix (compressed) to: {NORMALIZED_MATRIX_OUTPUT_TSV}")

    except Exception as e:
        print(f"An error occurred during export: {e}")


if __name__ == "__main__":
    export_data()