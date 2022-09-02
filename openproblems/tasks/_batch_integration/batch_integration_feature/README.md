<!--- TODO: add links --->

# Batch integration feature

This is a sub-task of the overall batch integration task. Batch (or data) integration
integrates datasets across batches that arise from various biological and technical
sources. Methods that integrate batches typically have three different types of output:
a corrected feature matrix, a joint embedding across batches, and/or an integrated
cell-cell similarity graph (e.g., a kNN graph). This sub-task focuses on all methods
that can output joint embeddings, and includes methods that canonically output corrected
feature matrices with subsequent postprocessing to generate a joint embedding. Other
sub-tasks for batch integration can be found for:

* [graphs](../batch_integration_graph/), and
* [embeddings](../batch_integration_embed/)

This sub-task was taken from a [benchmarking study of data integration
methods](https://www.biorxiv.org/content/10.1101/2020.05.22.111161v2).

## API

Datasets should contain the following attributes:

* `adata.obs["batch"]` with the batch covariate, and
* `adata.obs["label"]` with the cell identity label
* `adata.layers['counts']` with raw, integer UMI count data,
* `adata.layers['log_scran_pooling]` with log-normalized data and
* `adata.X` with log-normalized data

Methods should batch correct the matrix stored in `adata.X`.

The `openproblems-python-batch-integration` docker container is used for the methods
that
can be installed without package conflicts. For R methods, the `openproblems-r-extras`
container is used.

Methods are run in four different scenarios that include scaling and highly variable
gene selection:

* `full_unscaled`
* `hvg_unscaled`
* `full_scaled`
* `hvg_scaled`

Metrics can compare:

* `adata.X` vs `adata.layers['log_scran_pooling']`

To reuse metrics functions from `scIB`, [`metrics._utils._get_split`](metrics/_utils.py)
separates the combined anndata into an integrated and an unintegrated anndata object.


## Metrics
### HVG conservation
Metric that computes the average percentage of overlapping highly variable genes per
batch pre post integration.
