# Batch Integration with Embedding Output

This sub-task focuses on all methods that can output integrated embeddings.
Additionally, metrics can also be applied to PCA embeddings of the feature matrix output.
Other sub-tasks for batch integration can be found for:

* [graph](../graph/), and
* [corrected features](../feature/)

This sub-task was taken from
a [benchmarking study of data integration methods](https://www.biorxiv.org/content/10.1101/2020.05.22.111161v2).

## API

Datasets should contain the following attributes:

* `adata.uns['name']`: name of the dataset
* `adata.obs['batch']` with the batch covariate,
* `adata.obs['label']` with the cell identity label,
* `adata.layers['counts']` with raw, integer UMI count data, and
* `adata.X` with log-normalized data

Methods should assign output to:

* `adata.obsm['X_emb']`

Methods are run in four different scenarios that include scaling and highly variable gene selection:

* `full_unscaled`
* `hvg_unscaled`
* `full_scaled`
* `hvg_scaled`

Metrics can compare:

* `adata.obsm['X_emb']` to `adata.obsm['X_pca']`
* `adata.obsm['X_emb']` to `adata.obs['label']`
* `adata.obsm['X_emb']` to `adata.obs['batch']`
