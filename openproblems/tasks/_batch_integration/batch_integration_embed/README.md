<!--- TODO: add links --->

# Batch integration embedding

This is a sub-task of the overall batch integration task. Batch (or data) integration
integrates datasets across batches that arise from various biological and technical
sources. Methods that integrate batches typically have three different types of output:
a corrected feature matrix, a joint embedding across batches, and/or an integrated
cell-cell similarity graph (e.g., a kNN graph). This sub-task focuses on all methods
that can output joint embeddings, and includes methods that canonically output corrected
feature matrices with subsequent postprocessing to generate a joint embedding. Other
sub-tasks for batch integration can be found for:

* [graphs](../batch_integration_graph/), and
* [corrected features](../batch_integration_features)

This sub-task was taken from a
[benchmarking study of data integration
methods](https://www.biorxiv.org/content/10.1101/2020.05.22.111161v2).

## API

Datasets should contain the following attributes:

* `adata.obs["batch"]` with the batch covariate, and
* `adata.obs["label"]` with the cell identity label
* `adata.obsm['X_uni']` with a pre-integration embedding (usually PCA)
* `adata.layers['counts']` with raw, integer UMI count data, and
* `adata.X` with log-normalized data

Methods should assign output to `adata.obsm['X_emb'].

The `openproblems-python-batch-integration` docker container is used for the methods
that can be installed without package conflicts. For R methods, the
`openproblems-r-extras` container is used.

Most methods in this task are run in four different scenarios that include scaling and
highly variable gene selection:

* `full_unscaled`
* `hvg_unscaled`
* `full_scaled`
* `hvg_scaled`

Where `full` refers to the full gene set that is used as input to the method.

Metrics can compare:

* `adata.obsm['X_emb']` to `adata.obsm['X_uni']`,
* `adata.obsm['X_emb']` to `adata.obs['label']`, and/or
* `adata.obsm['X_emb']` to `adata.obs['batch']`.

To reuse metrics functions from `scIB`, [`metrics._utils._get_split`](metrics/_utils.py)
separates the combined anndata into an integrated and an unintegrated anndata object.
