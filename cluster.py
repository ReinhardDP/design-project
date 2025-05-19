import sys
import os
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import re
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
import spm1d

# Check command-line arguments
if len(sys.argv) != 3:
    print("Usage: python Clustering.py <adata_h5ad> <hotelling_csv>")
    sys.exit(1)

adata_path = sys.argv[1]
hotelling_path = sys.argv[2]

# pseudo-bulk expression matrix
adata = sc.read(adata_path, cache=False)
def compute_pseudo_bulk(adata, group_key='group'):
    """
    Calculate the pseudo-bulk expression matrix, i.e., the gene expression means grouped by `group_key`.
    """
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X  # Convert to dense matrices
    expr_df = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
    expr_df[group_key] = adata.obs[group_key].values
    pseudo_bulk = expr_df.groupby(group_key).mean()  # Aggregate gene expression by group
    return pseudo_bulk

pseudo_bulk = compute_pseudo_bulk(adata, group_key="group")
pseudo_bulk.to_csv("bulk_expression_matrix.csv")   # Save pseudo-bulk expression data
print("Pseudo-bulk expression matrix saved.")

# Correlation and clustering
bulk_corr = pseudo_bulk.T.corr(method='spearman')
dist_mat = pdist(bulk_corr, metric='cityblock')
linkage_matrix = linkage(dist_mat, method='complete')

# Load Hotelling T² results
hotelling = pd.read_csv(hotelling_path, index_col=0)
hotelling["score"] = pd.to_numeric(hotelling["HotellingT2"], errors='coerce')

# Get dendrogram order
ddata = dendrogram(linkage_matrix, no_plot=True, labels=bulk_corr.index.tolist())
original_order = ddata["ivl"]

# Sort variants by Hotelling T² score
sorted_variants = sorted(original_order, key=lambda v: hotelling["score"].get(v, 0), reverse=True)

# Cluster assignment
clusters = fcluster(linkage_matrix, t=3, criterion='maxclust')
variant_list = bulk_corr.index.tolist()
variant_to_cluster = {variant: cl for variant, cl in zip(variant_list, clusters)}

# Plot dendrogram
plt.figure(figsize=(20, 8))

dendrogram(
    linkage_matrix,
    labels=variant_list,
    leaf_rotation=90,
    leaf_font_size=12,
    color_threshold=0,
    above_threshold_color='black',
    link_color_func=lambda k: 'black'
)
plt.title("Hierarchical Clustering Dendrogram", fontsize=20)
plt.xticks(fontsize=10)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.savefig("word_dendrogram.png", dpi=300, bbox_inches='tight')
plt.savefig("latex_dendrogram.svg")
plt.close()
print("Hierarchical Clustering Dendrogram saved.")

# Save sorted variant order
# color mapping
color_map = {1: 'black', 2: 'lightblue', 3: 'blue'}
label_colors = {variant: color_map.get(variant_to_cluster.get(variant, 1), 'black') for variant in variant_list}
sorted_variants_df = pd.DataFrame({
    "variant": original_order,
    "color": [label_colors.get(v, 'black') for v in original_order]
})
sorted_variants_df.to_csv("sorted_variants.csv", index=False)
print("sorted_variants.csv saved.")
