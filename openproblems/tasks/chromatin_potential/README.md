# Chromatin accessibility prediction

Chromatin accessibility prediction refers to the gene expression prediction of a cell or cell type from ATAC-seq peaks. For a summary or all relevant models, see gene score method in [Jeffrey M. Granja et al.](https://www.biorxiv.org/content/10.1101/2020.04.28.066498v1), [Su Wang et al.](https://pubmed.ncbi.nlm.nih.gov/24263090/) et al.

## API

Datasets should contain the following attributes:

* `adata.obsm['mode2']` (cell by peak matrix of single cell atac-seq)
* `adata.X` (cell by gene matrix of single cell gene expression, which is the ground truth)
* `adata.obsm["gene_score"]` (model-based atac-seq prediction of gene regulation, which is the prediction)

Methods should assign gene regulation scores to `adata.obsm['gene_score']` using only single cell atac-seq peak counts. 

Metrics should compare adata.obsm['gene_score'] with the true gene expression in adata.X.

## Annotation
This task currently assumes Mus Musculus reference is used in the single cell dataset, which can be extended in the future if other species single cell sci-car, SNARE-seq or SHARE-seq is added.
