# Batch Integration with Graph Output

The output of all batch integration tasks can be represented as a graph. This sub-task focuses on all methods that can
output integrated graphs, and includes methods that canonically output the other two data formats with subsequent
postprocessing to generate a graph. Other sub-tasks for batch integration can be found for:

* [embeddings](../embedding/), and
* [corrected features](../feature/)

This sub-task was taken from
a [benchmarking study of data integration methods](https://www.biorxiv.org/content/10.1101/2020.05.22.111161v2).

## API

Datasets should contain the following attributes:

* `adata.obs['batch']` with the batch covariate,
* `adata.obs['label']` with the cell identity label,
* `adata.layers['counts']` with raw, integer UMI count data, and
* `adata.X` with log-normalized data

Methods should assign output to:

* `adata.obsp['connectivities']` and `adata.obsp['distances']`

Methods are run in four different scenarios that include scaling and highly variable gene selection:

* `full_unscaled`
* `hvg_unscaled`
* `full_scaled`
* `hvg_scaled`

Metrics can compare:

* `adata.obsp['connectivities']` to `adata.obs['uni_connectivies']`,
* `adata.obsp['connectivities']` to `adata.obs['label']`, and/or
* `adata.obsp['connectivities']` to `adata.obs['batch']`.
