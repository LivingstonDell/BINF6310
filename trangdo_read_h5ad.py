import scanpy as sc
import os

# The filename of the processed AnnData object
H5AD_FILE = "MBM_Christopher.h5ad"


def read_and_summarize():
    """
    Loads the processed AnnData object and prints a comprehensive summary
    of its contents, including dimensions and stored analysis results.
    """
    print("=" * 70)
    print(f"Attempting to read file: {H5AD_FILE}")
    print("=" * 70)

    if not os.path.exists(H5AD_FILE):
        print(f"Error: The input file '{H5AD_FILE}' was not found.")
        print("Please ensure the script is run in the same directory as the .h5ad file.")
        return

    try:
        # Load the AnnData object from the H5AD file
        adata = sc.read_h5ad(H5AD_FILE)

        print("\n--- Summary of Processed AnnData Object ---")
        # Printing the object automatically calls the summary function
        print(adata)

        print("\n--- Key Metadata Columns ---")
        # Show the first few columns of the cell metadata (obs)
        # This includes 'leiden', 'cell_types', 'pre_post', and 'outcome'
        print(adata.obs[['cell_types', 'outcome', 'pre_post', 'leiden']].head())

        print("\n--- Stored Dimensionality Reduction ---")
        # Show the keys for the calculated embeddings (UMAP, PCA)
        print("Available Embeddings (obsm keys):", adata.obsm_keys())


    except Exception as e:
        print(f"An error occurred while reading the H5AD file: {e}")


if __name__ == "__main__":
    read_and_summarize()
