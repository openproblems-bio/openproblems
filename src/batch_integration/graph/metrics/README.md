# Evaluation of Batch Integration with Graph Output

Metrics on graph output include:

* adjusted rand index ARI
* normalized mutual information NMI

## API

All datasets should contain the following attributes:

* `adata.uns['name']`: name of the dataset
* `adata.obs['batch']`: the batch covariate
* `adata.obs['label']`: the cell identity label
* `adata.obs['uni_connectivies']`: graph connectivities before integration
* `adata.obsp['connectivities']`: graph connectivities after integration
* `adata.obsp['distances']`: graph distances after integration

Metrics compare:

* `adata.obsp['connectivities']` to `adata.obs['uni_connectivies']`,
* `adata.obsp['connectivities']` to `adata.obs['label']`, and/or
* `adata.obsp['connectivities']` to `adata.obs['batch']`.
