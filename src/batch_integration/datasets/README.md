# Datasets

Viash component for preparing data **before** running data integration methods.

## API

This module creates Anndata objects that contain:

* `adata.uns['name']`: name of the dataset
* `adata.obs['batch']`: batch covariate
* `adata.obs['label']`: cell identity label
* `adata.vars['hvg']`: label whether a gene is identified as highly variable
* `adata.layers['counts']`: raw, integer UMI count data
* `adata.layers['logcounts']`: log-normalized count data
* `adata.layers['logcounts_scaled']`: scaled log-normalized count data
* `adata.X`: same as in `adata.layers['logcounts']`
