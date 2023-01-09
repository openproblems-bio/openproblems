<!--- TODO: add links --->

# Batch integration feature

## The task

This is a sub-task of the overall batch integration task. Batch (or data) integration
integrates datasets across batches that arise from various biological and technical
sources. Methods that integrate batches typically have three different types of output:
a corrected feature matrix, a joint embedding across batches, and/or an integrated
cell-cell similarity graph (e.g., a kNN graph). This sub-task focuses on all methods
that can output feature matrices. Other sub-tasks for batch integration can be found
for:

* [graphs](../batch_integration_graph/), and
* [embeddings](../batch_integration_embed/)

This sub-task was taken from a [benchmarking study of data integration
methods](https://openproblems.bio/bibliography#luecken2022benchmarking).

## The metrics

Metrics for batch integration (feature) measure how well feature-level information is
batch corrected. This is only done on by capturing biological variance conservation.
Further metrics for batch correction and biological variance conservation that are
calculated on lower dimensional feature spaces extrapolated from corrected feature
outputs can be found in the batch integration embed and graph tasks.

* **HVG conservation**: This metric computes the average percentage of overlapping
highly variable genes per batch before and after integration.

## API

WARNING: other than most tasks, `adata.X` should contain log-normalized data.
   This is the case as we are comparing the results of integration on the
   features pre- and post-integration and the data comes from different technologies.
   In this subtask, we are computing a pre-integration embedding on the normalized
   features.
   For UMI data, the data is scran-normalized, full-length data is TPM-normalized.

Datasets should contain the following attributes:

* `adata.obs["batch"]` with the batch covariate, and
* `adata.obs["label"]` with the cell identity label
* `adata.layers['counts']` with raw, integer UMI count data,
* `adata.layers['log_normalized']` with log-normalized data and
* `adata.X` with log-normalized data

Methods should store their a batch-corrected gene expression matrix in `adata.X`.

The `openproblems-python-batch-integration` docker container is used for the methods
that
can be installed without package conflicts. For R methods, the `openproblems-r-extras`
container is used.

Most methods in the current task are run in four different scenarios that include
caling
and highly variable gene selection:

* `full_unscaled`
* `hvg_unscaled`
* `full_scaled`
* `hvg_scaled`

Metrics for this task compare:

* `adata.X` (corrected adata) with `adata.layers['log_normalized']` (uncorrected data)

To reuse metrics functions from `scIB`, [`metrics._utils._get_split`](metrics/_utils.py)
separates the combined anndata into an integrated and an unintegrated anndata object.
