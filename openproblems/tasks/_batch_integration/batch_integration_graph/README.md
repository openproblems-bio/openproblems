<!--- TODO: add links --->

# Batch integration graph

This is a sub-task of the overall batch integration task. Batch (or data) integration methods integrate datasets across batches that arise from various biological and technical sources. Methods that integrate batches typically have three different types of output: a corrected feature matrix, a joint embedding across batches, and/or an integrated cell-cell similarity graph (e.g., a kNN graph). This sub-task focuses on all methods that can output integrated graphs, and includes methods that canonically output the other two data formats with subsequent postprocessing to generate a graph. Other sub-tasks for batch integration can be found for:

* [embeddings](../batch_integration_embed/), and
* [corrected features]()

This sub-task was taken from a [benchmarking study of data integration methods](https://www.biorxiv.org/content/10.1101/2020.05.22.111161v2).


## API

Datasets should contain the following attributes:

* `adata.obs["batch"]` with the batch covariate,
* `adata.obs["label"]` with the cell identity label,
* `adata.layers['counts']` with raw, integer UMI count data, and
* `adata.obsp['uni_connectivities']` with an unintegrated connectivity matrix generated
  by  `scanpy.pp.neighbors()`
* `adata.X` with log-normalized data

Methods can take anything from datasets as input and should assign output to:
* `adata.obsp['connectivities']` and `adata.obsp['distances']`, or
* `adata.uns['neighbors']['connectivities']` and  `adata.uns['neighbors']['distances']`, and

Please note, that most methods do not use cell type labels, which improves their usability.

The `openproblems-python-batch-integration` docker container is used for the methods that
can be installed without package conflicts. For R methods, the `openproblems-r-extras`
container is used.

Methods are run in four different scenarios that include scaling and highly variable gene selection:
* `full_unscaled`
* `hvg_unscaled`
* `full_scaled`
* `hvg_scaled`

An example script can be found [here](methods/_example.py)
Functions for scaling and highly variable gene selection per batch are reused from [`scIB`](https://github.com/theislab/scib). Additionally, method wrappers are reused from scIB where possible.

Metrics can compare:
* `adata.obsp['connectivities']` to `adata.obs['uni_connectivies']`,
* `adata.obsp['connectivities']` to `adata.obs['label']`, and/or
* `adata.obsp['connectivities']` to `adata.obs['batch']`.
