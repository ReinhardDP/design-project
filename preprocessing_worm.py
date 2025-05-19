import os
import sys
import pandas as pd
import anndata
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

# Check command-line arguments
if len(sys.argv) != 3:
    print("Usage: python qc_processing.py <input_matrix.csv> <metadata.csv>")
    sys.exit(1)

# Read input paths from arguments
input_matrix = sys.argv[1]
metadata_file = sys.argv[2]

# Hard-coded output parameters
figures_dir = "qc_figures"
output_h5ad = "GSE255865_processed.h5ad"
min_genes = 200
min_cells = 3

# Ensure the output directory exists
os.makedirs(figures_dir, exist_ok=True)

# 1. Load expression matrix (genes × samples)
expr = pd.read_csv(input_matrix, index_col=0)

# 2. Create AnnData object (transpose so samples become observations)
adata = anndata.AnnData(expr.T)
adata.var_names_make_unique()

# 3. Load and align metadata
meta = pd.read_csv(metadata_file).set_index("sampleID")
common = adata.obs_names.intersection(meta.index)
adata = adata[common]
adata.obs = meta.loc[common]

# 4. Compute QC metrics: total counts and number of genes detected per cell
adata.obs["n_counts"] = adata.X.sum(axis=1)
adata.obs["n_genes"]  = (adata.X > 0).sum(axis=1)

# 5. Plot and save QC distributions
plt.figure(figsize=(8, 5))
sns.histplot(adata.obs["n_counts"], bins=50, kde=True)
plt.title("Total Counts per Sample")
plt.savefig(os.path.join(figures_dir, "counts_per_sample.png"))
plt.close()

plt.figure(figsize=(8, 5))
sns.histplot(adata.obs["n_genes"], bins=50, kde=True)
plt.title("Detected Genes per Sample")
plt.savefig(os.path.join(figures_dir, "genes_per_sample.png"))
plt.close()

# 6. Filter cells and genes by QC thresholds
sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_genes(adata, min_cells=min_cells)

# 7. Save the processed AnnData to an h5ad file
adata.write_h5ad(output_h5ad)
print(f"Processed AnnData saved to {output_h5ad}")