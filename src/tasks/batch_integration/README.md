# Batch Integration

Remove unwanted batch effects from scRNA data while retaining
biologically meaningful variation.

Path:
[`src/tasks/batch_integration`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/tasks/batch_integration)

## Motivation

As single-cell technologies advance, single-cell datasets are growing
both in size and complexity. Especially in consortia such as the Human
Cell Atlas, individual studies combine data from multiple labs, each
sequencing multiple individuals possibly with different technologies.
This gives rise to complex batch effects in the data that must be
computationally removed to perform a joint analysis. These batch
integration methods must remove the batch effect while not removing
relevant biological information. Currently, over 200 tools exist that
aim to remove batch effects scRNA-seq datasets \[@zappia2018exploring\].
These methods balance the removal of batch effects with the conservation
of nuanced biological information in different ways. This abundance of
tools has complicated batch integration method choice, leading to
several benchmarks on this topic \[@luecken2020benchmarking;
@tran2020benchmark; @chazarragil2021flexible; @mereu2020benchmarking\].
Yet, benchmarks use different metrics, method implementations and
datasets. Here we build a living benchmarking task for batch integration
methods with the vision of improving the consistency of method
evaluation.

## Description

In this task we evaluate batch integration methods on their ability to
remove batch effects in the data while conserving variation attributed
to biological effects. As input, methods require either normalised or
unnormalised data with multiple batches and consistent cell type labels.
The batch integrated output can be a feature matrix, a low dimensional
embedding and/or a neighbourhood graph. The respective batch-integrated
representation is then evaluated using sets of metrics that capture how
well batch effects are removed and whether biological variance is
conserved. We have based this particular task on the latest, and most
extensive benchmark of single-cell data integration methods
\[@luecken2022benchmarking\].

## Authors & contributors

| name              | roles              |
|:------------------|:-------------------|
| Michaela Mueller  | maintainer, author |
| Kai Waldrant      | contributor        |
| Robrecht Cannoodt | contributor        |
| Daniel Strobl     | author             |

## API

``` mermaid
flowchart LR
  file_unintegrated(Unintegrated)
  file_integrated_embedding(Integrated embedding)
  file_integrated_feature(Integrated Feature)
  file_integrated_graaf(Integrated Graph)
  file_score(Score)
  file_common_dataset(Common dataset)
  comp_method_embedding[/Method (embedding)/]
  comp_method_feature[/Method (feature)/]
  comp_method_graaf[/Method (graph)/]
  comp_metric_embedding[/Metric (embedding)/]
  comp_metric_feature[/Metric (feature)/]
  comp_metric_graaf[/Metric (graph)/]
  comp_process_dataset[/Data processor/]
  file_unintegrated---comp_method_embedding
  file_unintegrated---comp_method_feature
  file_unintegrated---comp_method_graaf
  file_integrated_embedding---comp_metric_embedding
  file_integrated_feature---comp_metric_feature
  file_integrated_graaf---comp_metric_graaf
  file_common_dataset---comp_process_dataset
  comp_method_embedding-->file_integrated_embedding
  comp_method_feature-->file_integrated_feature
  comp_method_graaf-->file_integrated_graaf
  comp_metric_embedding-->file_score
  comp_metric_feature-->file_score
  comp_metric_graaf-->file_score
  comp_process_dataset-->file_unintegrated
```

## File format: Common dataset

A dataset processed by the common dataset processing pipeline.

Example file: `resources_test/common/pancreas/dataset.h5ad`

Description:

This dataset contains both raw counts and normalized data matrices, as
well as a PCA embedding, HVG selection and a kNN graph.

Format:

<div class="small">

    AnnData object
     obs: 'celltype', 'batch', 'tissue', 'size_factors'
     var: 'hvg', 'hvg_score'
     obsm: 'X_pca'
     obsp: 'knn_distances', 'knn_connectivities'
     varm: 'pca_loadings'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'dataset_name', 'data_url', 'data_reference', 'dataset_summary', 'dataset_description', 'dataset_organism', 'pca_variance', 'knn'

</div>

Slot description:

<div class="small">

| Slot                         | Type      | Description                                                                    |
|:-----------------------------|:----------|:-------------------------------------------------------------------------------|
| `obs["celltype"]`            | `string`  | (*Optional*) Cell type information.                                            |
| `obs["batch"]`               | `string`  | (*Optional*) Batch information.                                                |
| `obs["tissue"]`              | `string`  | (*Optional*) Tissue information.                                               |
| `obs["size_factors"]`        | `double`  | (*Optional*) The size factors created by the normalisation method, if any.     |
| `var["hvg"]`                 | `boolean` | Whether or not the feature is considered to be a ‘highly variable gene’.       |
| `var["hvg_score"]`           | `integer` | A ranking of the features by hvg.                                              |
| `obsm["X_pca"]`              | `double`  | The resulting PCA embedding.                                                   |
| `obsp["knn_distances"]`      | `double`  | K nearest neighbors distance matrix.                                           |
| `obsp["knn_connectivities"]` | `double`  | K nearest neighbors connectivities matrix.                                     |
| `varm["pca_loadings"]`       | `double`  | The PCA loadings matrix.                                                       |
| `layers["counts"]`           | `integer` | Raw counts.                                                                    |
| `layers["normalized"]`       | `double`  | Normalised expression values.                                                  |
| `uns["dataset_id"]`          | `string`  | A unique identifier for the dataset.                                           |
| `uns["dataset_name"]`        | `string`  | Nicely formatted name.                                                         |
| `uns["data_url"]`            | `string`  | (*Optional*) Link to the original source of the dataset.                       |
| `uns["data_reference"]`      | `string`  | (*Optional*) Bibtex reference of the paper in which the dataset was published. |
| `uns["dataset_summary"]`     | `string`  | Short description of the dataset.                                              |
| `uns["dataset_description"]` | `string`  | Long description of the dataset.                                               |
| `uns["dataset_organism"]`    | `string`  | (*Optional*) The organism of the sample in the dataset.                        |
| `uns["pca_variance"]`        | `double`  | The PCA variance objects.                                                      |
| `uns["knn"]`                 | `object`  | Supplementary K nearest neighbors data.                                        |

</div>

## Component type: Data processor

Path:
[`src/batch_integration`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/batch_integration)

A label projection dataset processor.

Arguments:

<div class="small">

| Name       | Type   | Description                                                    |
|:-----------|:-------|:---------------------------------------------------------------|
| `--input`  | `file` | A dataset processed by the common dataset processing pipeline. |
| `--output` | `file` | (*Output*) Unintegrated AnnData HDF5 file.                     |

</div>

## File format: Unintegrated

Unintegrated AnnData HDF5 file.

Example file:
`resources_test/batch_integration/pancreas/unintegrated.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     obs: 'batch', 'label'
     var: 'hvg'
     obsm: 'X_pca'
     obsp: 'knn_connectivities'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'normalization_id', 'dataset_organism'

</div>

Slot description:

<div class="small">

| Slot                         | Type      | Description                                                              |
|:-----------------------------|:----------|:-------------------------------------------------------------------------|
| `obs["batch"]`               | `string`  | Batch information.                                                       |
| `obs["label"]`               | `string`  | label information.                                                       |
| `var["hvg"]`                 | `boolean` | Whether or not the feature is considered to be a ‘highly variable gene’. |
| `obsm["X_pca"]`              | `double`  | The resulting PCA embedding.                                             |
| `obsp["knn_connectivities"]` | `double`  | K nearest neighbors connectivities matrix.                               |
| `layers["counts"]`           | `integer` | Raw counts.                                                              |
| `layers["normalized"]`       | `double`  | Normalized expression values.                                            |
| `uns["dataset_id"]`          | `string`  | A unique identifier for the dataset.                                     |
| `uns["normalization_id"]`    | `string`  | Which normalization was used.                                            |
| `uns["dataset_organism"]`    | `string`  | Which normalization was used.                                            |

</div>

## Component type: Method (embedding)

Path:
[`src/batch_integration/methods`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/batch_integration/methods)

A batch integration embedding method.

Arguments:

<div class="small">

| Name       | Type      | Description                                                                |
|:-----------|:----------|:---------------------------------------------------------------------------|
| `--input`  | `file`    | Unintegrated AnnData HDF5 file.                                            |
| `--output` | `file`    | (*Output*) An integrated AnnData HDF5 file.                                |
| `--hvg`    | `boolean` | (*Optional*) Whether to subset to highly variable genes. Default: `FALSE`. |

</div>

## Component type: Method (feature)

Path:
[`src/batch_integration/methods`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/batch_integration/methods)

A batch integration feature method.

Arguments:

<div class="small">

| Name       | Type      | Description                                                                |
|:-----------|:----------|:---------------------------------------------------------------------------|
| `--input`  | `file`    | Unintegrated AnnData HDF5 file.                                            |
| `--output` | `file`    | (*Output*) Integrated AnnData HDF5 file.                                   |
| `--hvg`    | `boolean` | (*Optional*) Whether to subset to highly variable genes. Default: `FALSE`. |

</div>

## Component type: Method (graph)

Path:
[`src/batch_integration/methods`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/batch_integration/methods)

A batch integration graph method.

Arguments:

<div class="small">

| Name       | Type      | Description                                                                |
|:-----------|:----------|:---------------------------------------------------------------------------|
| `--input`  | `file`    | Unintegrated AnnData HDF5 file.                                            |
| `--output` | `file`    | (*Output*) Integrated AnnData HDF5 file.                                   |
| `--hvg`    | `boolean` | (*Optional*) Whether to subset to highly variable genes. Default: `FALSE`. |

</div>

## File format: Integrated embedding

An integrated AnnData HDF5 file.

Example file: `resources_test/batch_integration/pancreas/scvi.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     obs: 'batch', 'label'
     var: 'hvg'
     obsm: 'X_pca', 'X_emb'
     obsp: 'knn_connectivities'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'normalization_id', 'dataset_organism', 'method_id', 'hvg', 'output_type'

</div>

Slot description:

<div class="small">

| Slot                         | Type      | Description                                                              |
|:-----------------------------|:----------|:-------------------------------------------------------------------------|
| `obs["batch"]`               | `string`  | Batch information.                                                       |
| `obs["label"]`               | `string`  | label information.                                                       |
| `var["hvg"]`                 | `boolean` | Whether or not the feature is considered to be a ‘highly variable gene’. |
| `obsm["X_pca"]`              | `double`  | The resulting PCA embedding.                                             |
| `obsm["X_emb"]`              | `double`  | integration embedding prediction.                                        |
| `obsp["knn_connectivities"]` | `double`  | K nearest neighbors connectivities matrix.                               |
| `layers["counts"]`           | `integer` | Raw counts.                                                              |
| `layers["normalized"]`       | `double`  | Normalized expression values.                                            |
| `uns["dataset_id"]`          | `string`  | A unique identifier for the dataset.                                     |
| `uns["normalization_id"]`    | `string`  | Which normalization was used.                                            |
| `uns["dataset_organism"]`    | `string`  | Which normalization was used.                                            |
| `uns["method_id"]`           | `string`  | A unique identifier for the method.                                      |
| `uns["hvg"]`                 | `boolean` | If the method was done on hvg or full.                                   |
| `uns["output_type"]`         | `string`  | what kind of output has been generated.                                  |

</div>

## File format: Integrated Feature

Integrated AnnData HDF5 file.

Example file: `resources_test/batch_integration/pancreas/combat.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     obs: 'batch', 'label'
     var: 'hvg'
     obsm: 'X_pca'
     obsp: 'knn_connectivities'
     layers: 'counts', 'normalized', 'corrected_counts'
     uns: 'dataset_id', 'normalization_id', 'dataset_organism', 'method_id', 'hvg', 'output_type'

</div>

Slot description:

<div class="small">

| Slot                         | Type      | Description                                                              |
|:-----------------------------|:----------|:-------------------------------------------------------------------------|
| `obs["batch"]`               | `string`  | Batch information.                                                       |
| `obs["label"]`               | `string`  | label information.                                                       |
| `var["hvg"]`                 | `boolean` | Whether or not the feature is considered to be a ‘highly variable gene’. |
| `obsm["X_pca"]`              | `double`  | The resulting PCA embedding.                                             |
| `obsp["knn_connectivities"]` | `double`  | K nearest neighbors connectivities matrix.                               |
| `layers["counts"]`           | `integer` | Raw counts.                                                              |
| `layers["normalized"]`       | `double`  | Normalized expression values.                                            |
| `layers["corrected_counts"]` | `double`  | Corrected counts after integration.                                      |
| `uns["dataset_id"]`          | `string`  | A unique identifier for the dataset.                                     |
| `uns["normalization_id"]`    | `string`  | Which normalization was used.                                            |
| `uns["dataset_organism"]`    | `string`  | Which normalization was used.                                            |
| `uns["method_id"]`           | `string`  | A unique identifier for the method.                                      |
| `uns["hvg"]`                 | `boolean` | If the method was done on hvg or full.                                   |
| `uns["output_type"]`         | `string`  | what kind of output has been generated.                                  |

</div>

## File format: Integrated Graph

Integrated AnnData HDF5 file.

Example file: `resources_test/batch_integration/pancreas/bbknn.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     obs: 'batch', 'label'
     var: 'hvg'
     obsm: 'X_pca'
     obsp: 'knn_connectivities', 'connectivities'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'normalization_id', 'dataset_organism', 'method_id', 'hvg', 'output_type'

</div>

Slot description:

<div class="small">

| Slot                         | Type      | Description                                                              |
|:-----------------------------|:----------|:-------------------------------------------------------------------------|
| `obs["batch"]`               | `string`  | Batch information.                                                       |
| `obs["label"]`               | `string`  | label information.                                                       |
| `var["hvg"]`                 | `boolean` | Whether or not the feature is considered to be a ‘highly variable gene’. |
| `obsm["X_pca"]`              | `double`  | The resulting PCA embedding.                                             |
| `obsp["knn_connectivities"]` | `double`  | K nearest neighbors connectivities matrix.                               |
| `obsp["connectivities"]`     | `double`  | Neighbors connectivities matrix.                                         |
| `layers["counts"]`           | `integer` | Raw counts.                                                              |
| `layers["normalized"]`       | `double`  | Normalized expression values.                                            |
| `uns["dataset_id"]`          | `string`  | A unique identifier for the dataset.                                     |
| `uns["normalization_id"]`    | `string`  | Which normalization was used.                                            |
| `uns["dataset_organism"]`    | `string`  | Which normalization was used.                                            |
| `uns["method_id"]`           | `string`  | A unique identifier for the method.                                      |
| `uns["hvg"]`                 | `boolean` | If the method was done on hvg or full.                                   |
| `uns["output_type"]`         | `string`  | what kind of output has been generated.                                  |

</div>

## Component type: Metric (embedding)

Path:
[`src/batch_integration/metrics`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/batch_integration/metrics)

A batch integration embedding metric.

Arguments:

<div class="small">

| Name                 | Type   | Description                      |
|:---------------------|:-------|:---------------------------------|
| `--input_integrated` | `file` | An integrated AnnData HDF5 file. |
| `--output`           | `file` | (*Output*) Metric score file.    |

</div>

## Component type: Metric (feature)

Path:
[`src/batch_integration/metrics`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/batch_integration/metrics)

A batch integration feature metric.

Arguments:

<div class="small">

| Name                 | Type   | Description                   |
|:---------------------|:-------|:------------------------------|
| `--input_integrated` | `file` | Integrated AnnData HDF5 file. |
| `--output`           | `file` | (*Output*) Metric score file. |

</div>

## Component type: Metric (graph)

Path:
[`src/batch_integration/metrics`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/batch_integration/metrics)

A batch integration graph metric.

Arguments:

<div class="small">

| Name                 | Type   | Description                   |
|:---------------------|:-------|:------------------------------|
| `--input_integrated` | `file` | Integrated AnnData HDF5 file. |
| `--output`           | `file` | (*Output*) Metric score file. |

</div>

## File format: Score

Metric score file

Example file: `score.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     uns: 'dataset_id', 'normalization_id', 'method_id', 'metric_ids', 'metric_values', 'hvg', 'output_type'

</div>

Slot description:

<div class="small">

| Slot                      | Type      | Description                                                                                  |
|:--------------------------|:----------|:---------------------------------------------------------------------------------------------|
| `uns["dataset_id"]`       | `string`  | A unique identifier for the dataset.                                                         |
| `uns["normalization_id"]` | `string`  | Which normalization was used.                                                                |
| `uns["method_id"]`        | `string`  | A unique identifier for the method.                                                          |
| `uns["metric_ids"]`       | `string`  | One or more unique metric identifiers.                                                       |
| `uns["metric_values"]`    | `double`  | The metric values obtained for the given prediction. Must be of same length as ‘metric_ids’. |
| `uns["hvg"]`              | `boolean` | If the method was done on hvg or full.                                                       |
| `uns["output_type"]`      | `string`  | what kind of output has been generated.                                                      |

</div>
