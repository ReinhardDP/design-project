import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import sys
import seaborn as sns

# Argument parsing
if len(sys.argv) != 4:
    print("Usage: python script.py <input_file> <figures_dir> <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
figures_dir = sys.argv[2]
output_file = sys.argv[3]
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80, color_map='viridis')
sc.logging.print_versions()

# Load the data
adata = sc.read(input_file, cache=False)

# Make a box plot to see data before normalization
# Convert expression matrix to DataFrame
# If adata.X is sparse, convert to dense first
expression_df = pd.DataFrame(adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X,
                             index=adata.obs_names,
                             columns=adata.var_names)

subset_df = expression_df.iloc[:, :100]
# Melt the DataFrame for seaborn
melted_df = subset_df.melt(var_name='Gene', value_name='Log(Expression + 1)')

# BOX Plot
plt.figure(figsize=(12, 6))
sns.boxplot(x='Gene', y='Log(Expression + 1)', hue='Gene', data=melted_df, showfliers=False, palette="Set2", dodge=False)
plt.legend([],[], frameon=False)  # Hide legend
plt.title('Box Plot of the first 100 Log-Normalized Gene Expression Levels')
plt.xticks(rotation=90, fontsize=7)
plt.savefig(f'{figures_dir}/boxplot_pre.png',bbox_inches = "tight")
plt.close()

# Normalize
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
adata.raw = sc.pp.log1p(adata, copy=True)
sc.pp.log1p(adata)

# Select highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=4, min_disp=0.5, flavor='seurat')
adata = adata[:, adata.var['highly_variable']]

# Scale and PCA
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
adata.obsm['X_pca'] *= -1

variance = adata.uns['pca']['variance']
total_variance = sum(variance)
variance = (variance/total_variance)*100
cumulative_variance = np.cumsum(variance)

# PCA Plots
plt.figure(figsize=(8,5))
plt.plot(range(1, len(cumulative_variance)+1), cumulative_variance, marker='o', linestyle='--', color='b')
plt.xlabel('Number of Principal Components')
plt.ylabel('Cumulative Variance Explained in %')
plt.title('Cumulative Variance Explained by PCA')
plt.grid(True)
plt.savefig(f'{figures_dir}/pca_cumulative_variance.png', dpi=300, bbox_inches='tight')
plt.close()


plt.figure(figsize=(8,5))
plt.plot(range(1, len(variance)+1), variance, marker='o', linestyle='--', color='b')
plt.xlabel('Number of Principal Components')
plt.ylabel('Variance Explained in %')
plt.title('Variance Explained by PCA')
plt.grid(True)
plt.savefig(f'{figures_dir}/pca_variance.png', dpi=300, bbox_inches='tight')
plt.close()

# Make a box plot to see data after normalization
expression_df = pd.DataFrame(adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X,
                             index=adata.obs_names,
                             columns=adata.var_names)

subset_df = expression_df.iloc[:, :100]
# Melt the DataFrame for seaborn
melted_df = subset_df.melt(var_name='Gene', value_name='Log(Expression + 1)')

# BOX Plot
plt.figure(figsize=(12, 6))
sns.boxplot(x='Gene', y='Log(Expression + 1)', hue='Gene', data=melted_df, showfliers=False, palette="Set2", dodge=False)
plt.legend([],[], frameon=False)
plt.title('Box Plot of the first 100 Log-Normalized Gene Expression Levels')
plt.xticks(rotation=90, fontsize=7)
plt.savefig(f'{figures_dir}/boxplot_post.png',bbox_inches = "tight")
plt.close()

# Save processed data
adata.write(output_file)

print(f"Processed data saved to {output_file}")
