# Preprocessing_TP53
## Usage
```
python script.py <input_matrix> <input_genes> <input_cells> <variant2cell> <figures_dir> <output_file>
```

## Input Files

### 1. `input_matrix` (Matrix Market file)
- **Format:** `.mtx`
- **Description:** A sparse matrix in Matrix Market format representing gene expression counts.
- **Structure:**
  - Rows correspond to genes.
  - Columns correspond to cells.
  - Entries represent the expression levels.
- **Example:**
  ```
  %%MatrixMarket matrix coordinate integer general
  3 3 4
  1 1 5
  1 2 3
  2 2 8
  3 3 6
  ```

### 2. `input_genes` (Gene names file)
- **Format:** `.txt` or `.csv`
- **Description:** A file containing gene names, one per line, corresponding to the rows of `input_matrix`.
- **Structure:**
  - Single-column text file.
- **Example:**
  ```
  GeneA
  GeneB
  GeneC
  ```

### 3. `input_cells` (Cell names file)
- **Format:** `.txt` or `.csv`
- **Description:** A file containing cell identifiers, one per line, corresponding to the columns of `input_matrix`.
- **Structure:**
  - Single-column text file.
- **Example:**
  ```
  Cell1
  Cell2
  Cell3
  ```

### 4. `variant2cell` (Variant-to-cell mapping file)
- **Format:** `.csv` or `.txt`
- **Description:** A mapping file that links genetic variants to specific cells.
- **Structure:**
  - Two-column format with a header.
  - First column: Variant ID.
  - Second column: Cell identifier.
- **Example:**
  ```
  VariantID,CellID
  rs123,Cell1
  rs456,Cell2
  ```
### 5. `figures_dir` 
- **Description:** The directory where the figures will be stored.

## Output File
### `output_file`
- **Format:** `.h5ad`
- **Description:** Processed AnnData object saved in HDF5 format.

# Normalize

## Overview
This script processes single-cell RNA sequencing (scRNA-seq) data by performing normalization, variance filtering, principal component analysis (PCA), and visualization. It outputs processed data in AnnData format (`.h5ad`) and generates figures for quality control.

## Usage
```
python script.py <input_file> <figures_dir> <output_file>
```

### Arguments
- `<input_file>`: Path to the input AnnData (`.h5ad`) file containing raw scRNA-seq data.
- `<figures_dir>`: Path where the figures will be saved.
- `<output_file>`: Path where the processed AnnData (`.h5ad`) file will be saved.

## Input File
- The input file should be in `.h5ad` format, a standard format for storing single-cell RNA sequencing data.
- It should contain:
  - `adata.X`: Gene expression matrix (cells × genes)
  - `adata.obs_names`: Cell identifiers
  - `adata.var_names`: Gene names

## Processing Steps
1. **Load the Data**
   - Reads the input `.h5ad` file using `scanpy`.
   
2. **Pre-Normalization Quality Control**
   - Generates a boxplot (`figures/boxplot_pre.png`) for the first 100 genes to visualize expression levels before normalization.

3. **Normalization and Log Transformation**
   - Normalizes expression per cell.
   - Applies log1p transformation.
   
4. **Highly Variable Gene Selection**
   - Filters genes based on variance to retain only highly variable genes.

5. **Scaling and Principal Component Analysis (PCA)**
   - Scales gene expression values.
   - Performs PCA with 50 components.
   - Saves PCA variance and cumulative variance plots:
     - `figures/pca_variance.png`
     - `figures/pca_cumulative_variance.png`
   
6. **Post-Normalization Quality Control**
   - Generates a boxplot (`figures/boxplot_post.png`) for the first 100 genes to visualize expression levels after normalization.

7. **Save Processed Data**
   - Writes the processed data to the specified output `.h5ad` file.

## Dependencies
Ensure the following Python packages are installed:
```
pip install numpy pandas scanpy matplotlib seaborn igraph louvain
```

## Output Files
- Processed `.h5ad` file.
- Boxplots before and after normalization (`figures/boxplot_pre.png`, `figures/boxplot_post.png`).
- PCA variance plots (`figures/pca_variance.png`, `figures/pca_cumulative_variance.png`).

# Scoring
## Usage
```
python script.py <input_path> <reference_variant> <figures_directory> <output_file>
```
## Input arguments
### 1. `input path` 
- Path to the preprocessed and normalised .h5ad file

### 2. `reference variant`
- The reference variant against which all other variants will be compared. This argument has to be a (part of a) string that occurs in the name of each reference sample.
- Example: ctrl
  
### 3. `figures directory`
- Directory where the output figure will be saved.
  
## Output file
- CSV file that will store the Hotelling's T² results.
- Example output:
```
  variant,result
  VariantA,3.45
  VariantB,5.67
```

# Clustering
## Usage
```
python Cluster.py <adata_processed.h5ad> <HotellingT2.csv>
```
## Input arguments
- Path to the adata_processed.h5ad and HotellingT2.csv
  
## Output files
- Four clustering plots, bulk_expression_matrix.csv and TP53.corrL1.sorted_variants.csv  
- Example output:
  bulk_expression_matrix.csv:
  ```
  Rows: Variant groups (values of the variant column)
  Columns: Genes (adata.var_names)
  Values: Mean expression level of each gene across all cells in the corresponding variant group
  ```
  TP53.corrL1.sorted_variants.csv:
  ```
  Variant1 Black
  Variant2 Black
  Variant3 Blue
  ```
