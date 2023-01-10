<!--- TODO: add links --->

# Batch integration (graph)

## The task

This is a sub-task of the overall batch integration task. Batch (or data) integration
methods integrate datasets across batches that arise from various biological and
technical sources. Methods that integrate batches typically have three different types
of output: a corrected feature matrix, a joint embedding across batches, and/or an
integrated cell-cell similarity graph (e.g., a kNN graph). This sub-task focuses on all
methods that can output integrated graphs, and includes methods that canonically output
the other two data formats with subsequent postprocessing to generate a graph. Other
sub-tasks for batch integration can be found for:

* [embeddings](../batch_integration_embed/), and
* [corrected features](../batch_integration_feature/)

This sub-task was taken from a [benchmarking study of data integration
methods](https://openproblems.bio/bibliography#luecken2022benchmarking).

## The metrics

Metrics for batch integration (graph) measure how well batches are mixed while
biological signals are preserved. They are divided into batch correction and biological
variance conservation metrics.

### Batch correction

* **Graph connectivity**: The graph connectivity metric assesses whether the kNN graph
representation, G, of the integrated data connects all cells with the same cell identity
label.

### Biological variance removal

* **Adjusted rand index (ARI)**: The Rand index compares the overlap of two clusterings;
it considers both correct clustering overlaps while also counting correct disagreements
between two clusterings.
* **Iso label F1 score**: Isolated cell labels are identified as the labels present in
the least number of batches in the integration task. The score evaluates how well these
isolated labels separate from other cell identities based on clustering.
* **Normalized mutual information (NMI)**: NMI compares the overlap of two clusterings.
We used NMI to compare the cell-type labels with Louvain clusters computed on the
integrated dataset.

## API

WARNING: other than most tasks, `adata.X` should contain log-normalized data.
   This is the case as we are comparing the results of integration on the normalized
   features pre- and post-integration and the data comes from different technologies.
   In this subtask, we are computing a pre-integration graph on the normalized
   features.
   For UMI data, the data is scran-normalized, full-length data is TPM-normalized.

Datasets should contain the following attributes:

* `adata.obs["batch"]` with the batch covariate,
* `adata.obs["label"]` with the cell identity label,
* `adata.layers['counts']` with raw, integer UMI count data, and
* `adata.obsm['X_uni']` with the PCA embedding of the unintegrated representation
* `adata.obsp['uni_connectivities']` with an unintegrated connectivity matrix generated
  by  `scanpy.pp.neighbors()`
* `adata.X` with log-normalized data
* `adata.layers['log_normalized']` with log-normalized data

Methods can take anything from datasets as input and should assign output to:

* `adata.obsp['connectivities']` and `adata.obsp['distances']`, or
* `adata.uns['neighbors']['connectivities']` and  `adata.uns['neighbors']['distances']`.

Please note, that most methods do not use cell type labels, which improves their
usability.

The `openproblems-python-batch-integration` docker container is used for the methods
that can be installed without package conflicts. (NOTE: add additional containers here)
For R methods, the `openproblems-r-extras` container is used.

Methods are run in four different scenarios that include scaling and highly variable
gene selection:

* `full_unscaled`
* `hvg_unscaled`
* `full_scaled`
* `hvg_scaled`

Functions for scaling and highly variable gene selection per batch are reused from
[`scib`](https://github.com/theislab/scib). Additionally, method wrappers are reused
from scIB where possible.

Metrics can compare:

* `adata.obsp['connectivities']` to `adata.obs['uni_connectivies']`,
* `adata.obsp['connectivities']` to `adata.obs['label']`, and/or
* `adata.obsp['connectivities']` to `adata.obs['batch']`.
