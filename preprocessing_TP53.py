import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import sys
import seaborn as sns

# Argument parsing
if len(sys.argv) != 7:
    print("Usage: python script.py <input_matrix> <input_genes> <input_cells> <variant2cell> <figures_dir> <output_file>")
    sys.exit(1)

input_matrix = sys.argv[1]
input_genes = sys.argv[2]
input_cells = sys.argv[3]
variant2cell_file = sys.argv[4]
figures_dir = sys.argv[5]
output_file = sys.argv[6]

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80, color_map='viridis')
sc.logging.print_versions()

# Load the data
adata = sc.read(input_matrix, cache=True)
adata.var_names = pd.read_csv(input_genes, header=None)[0]
adata.obs_names = pd.read_csv(input_cells, header=None)[0]
adata.var_names_make_unique()

# Load the variant2cell.csv file
variant2cell = pd.read_csv(variant2cell_file, delimiter='\t')

variant2cell['variant'] = variant2cell['variant.detailed_multi']

# Set the 'cell' column as the index to align with adata.obs_names
variant2cell = variant2cell.set_index('cell')

# Ensure the indices in variant2cell match adata.obs_names
common_indices = variant2cell.index.intersection(adata.obs_names)
adata = adata[common_indices]  # Subset adata to match variant2cell
variant2cell = variant2cell.loc[common_indices]  # Subset variant2cell to match adata

# Add the metadata to adata.obs
adata.obs = pd.concat([adata.obs, variant2cell], axis=1)

# Remove cells with less than 200 genes and genes present in less than 3 cells
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)


# Add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

# Plot histogram of total counts per cell
plt.figure(figsize=(8, 5))
sns.histplot(adata.obs['n_counts'], bins=50, kde=True, color='skyblue')
plt.xlabel("Total Counts per Cell")
plt.ylabel("Number of Cells")
plt.title("Distribution of Total Counts per Cell")
plt.savefig(f'{figures_dir}/counts_per_cell.png',bbox_inches = "tight")
plt.close()

# Detected genes per cell (non-zero counts)
adata.obs['n_genes'] = (adata.X > 0).sum(axis=1).A1

# Plot histogram of total genes per cell
plt.figure(figsize=(8, 5))
sns.histplot(adata.obs['n_genes'], bins=50, kde=True, color='lightgreen')
plt.xlabel("Number of Detected Genes per Cell")
plt.ylabel("Number of Cells")
plt.title("Distribution of Detected Genes per Cell")
plt.savefig(f'{figures_dir}/genes_per_cell.png',bbox_inches = "tight")
plt.close()


# Save processed data
adata.write(output_file)

print(f"Processed data saved to {output_file}")
