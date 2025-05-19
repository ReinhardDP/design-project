import scanpy as sc
import pandas as pd
import spm1d
import sys
import matplotlib.pyplot as plt
import numpy as np
import csv

if len(sys.argv) != 5:
    print("Usage: python script.py <input_path> <reference_variant> <figures_directory> <output_file>")
    sys.exit(1)

input_path = sys.argv[1]
reference = sys.argv[2]
figures_dir = sys.argv[3]
output_file = sys.argv[4]

# Load AnnData object
adata = sc.read(input_path, cache=False)

#adata2 = sc.read('/data/dp6/lobvdrsc/adata_processed.h5ad', cache=False)
#cells_to_keep =set(adata2.obs_names)

variance = adata.uns['pca']['variance']
total_variance = sum(variance)
variance = (variance/total_variance)*100

num_pcs = 0
cumulative_variance = 0
for i in variance:
    cumulative_variance += i
    num_pcs += 1
    if cumulative_variance >= 80:
        break

print(f"{num_pcs} principal components were selected with a total explained variance of {cumulative_variance}%")

# Extract variant and PCA columns
old_variants = [cat for cat in adata.obs["variant"].cat.categories if reference not in cat]
print('old: ', old_variants) 
variants = [
    word.strip() 
    for cat in adata.obs["variant"].cat.categories 
    if reference not in cat
    for word in cat.split(',')
]
print('new: ', variants)
#pcas = adata2.obs.columns.tolist()[121:142]
adata_pca = adata.obsm["X_pca"][0:num_pcs]
#reference = 'P359P'

pc_names = [f"PC{i}" for i in range(num_pcs)]
adata_pca=pd.DataFrame(adata.obsm['X_pca'][:,:num_pcs],
                     index=adata.obs_names,
                    columns=pc_names)



# Filter df (PCA dataset)
#adata_pca = adata_pca[adata_pca.index.isin(cells_to_keep)].copy()

# Filter adata (AnnData object)
#adata = adata[adata.obs_names.isin(cells_to_keep)].copy()

results = []
score_values = []

referencedata = adata_pca[adata.obs['variant'].str.contains(reference, case=False, na=False, regex=False)]

for group in variants:
    variantdata = adata_pca[adata.obs['variant'].str.contains(group, case=False, na=False, regex=False)]

    T2 = spm1d.stats.hotellings2(referencedata, variantdata)
    value = T2.z
    results.append((group, value))
    score_values.append(value)



# Write results to a CSV file
with open(output_file, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["variant", "result"])  
    for variant, result in results:
        writer.writerow([variant, result]) 

# Plot histogram
plt.figure(figsize=(10, 6))
plt.hist(score_values, bins=20, color='skyblue', edgecolor='black')
plt.title("Histogram of Hotelling's T² Scores")
plt.xlabel("T² Score (z)")
plt.ylabel("Frequency")
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig(f"{figures_dir}/scores_histogram.png")
plt.close()
