# Chromatin accessibility prediction

Chromatin accessibility prediction refers to the gene expression prediction of a cell or
cell type from ATAC-seq peaks. For a summary or all relevant models, see gene score
method in [Jeffrey M. Granja et
al.](https://www.biorxiv.org/content/10.1101/2020.04.28.066498v1), [Su Wang et
al.](https://pubmed.ncbi.nlm.nih.gov/24263090/) et al.

## API

Datasets should contain the following attributes:

* `adata.uns['species']` (ensembl species name, e.g. `"mus_musculus"`)
* `adata.uns['release']` (ensembl release, e.g. `"100"`)
* `adata.uns['mode2_var_chr']` (single cell atac-seq peak chromosome)
* `adata.uns['mode2_var_start']` (single cell atac-seq peak start position)
* `adata.uns['mode2_var_end']` (single cell atac-seq peak end position)
* `adata.obsm['mode2']` (cell by peak matrix of single cell atac-seq)
* `adata.X` (cell by gene matrix of single cell gene expression, which is the ground
  truth)

Methods should assign gene regulation scores to `adata.obsm['gene_score']` using only
single cell atac-seq peak counts.

Metrics should compare `adata.obsm['gene_score']` with the true gene expression in
`adata.X`.
