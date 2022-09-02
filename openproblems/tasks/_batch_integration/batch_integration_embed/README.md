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

Methods should assign output to `adata.obsm['X_emb']`.

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

## Metrics

In the following, we will give a short description of the implemented metrics. We split by metrics capturing batch correction meaning the removal of batch effects and metrics describing biological conservation, meaning how well the biological differences between cell states are conserved.

### Batch correction metrics

#### kBET

The kBET algorithm (v.0.99.6, release 4c9dafa) determines whether the label composition
of a k nearest neighborhood of a cell is similar to the expected (global) label
composition (Buettner et al., Nat Meth 2019). The test is repeated for a random subset
of cells, and the results are summarized as a rejection rate over all tested
neighborhoods.

#### Silhouette batch score

We consider the absolute silhouette width, s(i), on
batch labels per cell i. Here, 0 indicates that batches are well mixed, and any
deviation from 0 indicates a batch effect.

#### Principal component regression

Compare the explained variance before and after integration. Return  a score between 0 and 1 (scaled=True) with 0 if the variance contribution hasn’t changed. The larger the score, the more different the variance contributions are before and after integration.

### Biological conservation metrics

#### Cell cycle score

The cell-cycle conservation score evaluates how well the cell-cycle effect can be
captured before and after integration.

#### Isolated label silhouette

This score evaluates for each cell type label that is present in more than one batch how isolated the cells of said label are from other cells.

#### Cell type ASW

For the bio-conservation score, the ASW was computed on cell identity labels.
