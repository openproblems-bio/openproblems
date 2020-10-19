# Chromatin accessibility prediction

Chromatin accessibility prediction refers to the gene expression prediction of a cell or cell type from ATAC-seq peaks. For a summary or all relevant models, see gene score method in [Jeffrey M. Granja et al.](https://www.biorxiv.org/content/10.1101/2020.04.28.066498v1), [Su Wang et al.](https://pubmed.ncbi.nlm.nih.gov/24263090/) et al.

## API

Datasets should contain the following attributes:

* `adata.obsm["gene_score"]` (model-based atac-seq prediction of gene regulation)
* `adata.obs["atac_rna_cor"]` (predicted correlation between gene expression and atac-seq)

Methods should assign gene regulation scores to `adata.obs['gene_score']` using only single cell atac-seq peak counts. 

Metrics should maximize the median correlation metric stored in `adata['atac_rna_cor']`.
