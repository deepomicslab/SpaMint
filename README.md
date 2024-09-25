# SpaMint
The software implementation of the method in 
[Deciphering more accurate cell-cell interactions by modeling cells and their interactions]().

<img src="https://github.com/deepomicslab/SpexMod/raw/main/main.pg" alt="Spint-Main">

# Pre-requirements
* numpy, pandas==1.5.2
* scipy, scanpy, umap
* loess
* smurf-imputation
  
# Installation
<details><summary>Expand section</summary>
  
## Installation with pip
```shell
pip install pySPAMINT
```
## Installation from the source code
conda create --name spamint python=3.9
conda activate spamint
```shell
wget https://github.com/deepomicslab/SpaMint/archive/refs/heads/main.zip
unzip main.zip
cd SpaMint-main
python setup.py install
```
 </details>


# Input file format
## 1. DataFrame format
<details><summary>Expand section</summary>

* **Spatial Transcriptomics (ST) Count Data**
  * `st_exp` dataframe with spots as rows and genes as columns
 
* **Spatial coordinates**
  * `st_coord` dataframe with spot as rows, axis x and y as columns 

* **Cell-type deconvoluted spatial matrix**
  * `st_decon` dataframe with spot as rows and cell-type as columns

* **Single-cell RNA-seq Count Data**
  * `sc_exp` dataframe with cells as rows and genes as columns

* **Single-cell RNA-seq Metadata**
  * `sc_meta` dataframe with cells as rows and cell types as columns
  * `cell_type_key` column name of the celltype identity in `sc_meta`

* **Single-cell RNA-seq distribution Data**
  * `sc_distribution` dataframe with cells as rows and genes as columns
    
* **Ligand and Receptor Data (optional)**
  * `lr_df` user provided dataframe with ligand-receptor pairs as rows, ligand, receptor and its weight as columns

***
Convert to adata format
```python
sc_adata, st_adata, sc_distribution, lr_df = pp.prep_adata(sc_exp = sc_exp, st_exp = st_exp, sc_distribution = sc_smurf, 
                            sc_meta = sc_meta, st_coord = st_coord, SP = species)
```
</details>

## 2. Adata format
<details><summary>Expand section</summary>
  
* **Spatial Transcriptomics (ST) Count Data**
  * `st_adata` adata.X with spots as rows and genes as columns
  * `st_adata.obs`  dataframe with spot as rows, spot coordinates x and y as columns 
 
* **Cell-type deconvoluted spatial matrix**
  * `st_decon` dataframe with spot as rows and cell-type as columns

* **Single-cell RNA-seq Count Data**
  * `sc_adata` adata.X dataframe with cells as rows and genes as columns
  * `sc_adata.obs` dataframe with cells as rows and cell types as columns

* **Single-cell RNA-seq distribution Data**
  * `sc_distribution` dataframe with cells as rows and genes as columns
</details>

# Usages
## Prep object
<details><summary>Expand section</summary>

```python
obj = spamint.spaMint(save_path = outDir, st_adata = st_adata, weight = st_decon, 
                 sc_distribution = sc_distribution, sc_adata = sc_adata, cell_type_key = 'celltype', 
                 st_tp = st_tp)
obj.prep()
```
### Parameters
* `save_path` Output Dir to save results
  
* `st_adata` adata.X Spatial Transcriptomics (ST) Count Data with spots as rows and genes as columns
  * `st_adata.obs`  dataframe with spot as rows, spot coordinates x and y as columns
    
* `weight` Cell-type deconvoluted spatial dataframe with spot as rows and cell-type as columns
    
* `sc_distribution` Single-cell RNA-seq distribution dataframe with cells as rows and genes as columns
    
* `sc_adata` adata.X Single-cell RNA-seq Count dataframe with cells as rows and genes as columns
  * `sc_adata.obs` dataframe with cells as rows and cell types as columns
 
* `cell_type_key` cell type colname in sc_adata.obs

* `st_tp` ST sequencing platform choose from st (ST legacy), visium (10X Visium), or slide-seq (Any single-cell resolution data)


</details>

## Cell selection
<details><summary>Expand section</summary>
select_cells(self, p = 0.1, mean_num_per_spot = 10,  max_rep = 3, repeat_penalty = 10)
        - p: Persentage of interface similarity for cell selection.
        - mean_num_per_spot: Average number of cells per spot.
        - max_rep: Maximum number of repetitions for cell selection.
        - repeat_penalty: Penalty applied for repeated selections.

gradient_descent(self, alpha, beta, gamma, delta, eta, 
                init_sc_embed = False,
                iteration = 20, k = 2, W_HVG = 2,
                left_range = 1, right_range = 2, steps = 1, dim = 2)

        Parameters:
        - alpha, beta, gamma, delta: Hyperparameters for the loss function.
        - eta: Learning rate for gradient descent.
        - init_sc_embed: Initial embedding for single-cell data.
        - iteration: Number of iterations for optimization.
        - k: Number of neighbors for KNN calculations.
        - W_HVG: Weight for highly variable genes.
        - left_range, right_range: Range parameters for embeddings.
        - steps: Number of steps for embedding adjustments.
        - dim: Dimensionality for embedding.
</details>
