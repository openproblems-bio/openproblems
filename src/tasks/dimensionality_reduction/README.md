# Dimensionality reduction for visualization

Reduction of high-dimensional datasets to 2D for visualization &
interpretation

Path:
[`src/tasks/dimensionality_reduction`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/tasks/dimensionality_reduction)

## Motivation

Dimensionality reduction is one of the key challenges in single-cell
data representation. Routine single-cell RNA sequencing (scRNA-seq)
experiments measure cells in roughly 20,000-30,000 dimensions (i.e.,
features - mostly gene transcripts but also other functional elements
encoded in mRNA such as lncRNAs). Since its inception,scRNA-seq
experiments have been growing in terms of the number of cells measured.
Originally, cutting-edge SmartSeq experiments would yield a few hundred
cells, at best. Now, it is not uncommon to see experiments that yield
over [100,000 cells](https://www.nature.com/articles/s41586-018-0590-4)
or even [\> 1 million cells](https://doi.org/10.1126/science.aba7721).

## Description

Each *feature* in a dataset functions as a single dimension. While each
of the ~30,000 dimensions measured in each cell contribute to an
underlying data structure, the overall structure of the data is
challenging to display in few dimensions due to data sparsity and the
[*“curse of
dimensionality”*](https://en.wikipedia.org/wiki/Curse_of_dimensionality)
(distances in high dimensional data don’t distinguish data points well).
Thus, we need to find a way to [dimensionally
reduce](https://en.wikipedia.org/wiki/Dimensionality_reduction) the data
for visualization and interpretation.

## Authors & contributors

| name                   | roles              |
|:-----------------------|:-------------------|
| Luke Zappia            | maintainer, author |
| Michal Klein           | author             |
| Scott Gigante          | author             |
| Ben DeMeo              | author             |
| Juan A. Cordero Varela | contributor        |
| Robrecht Cannoodt      | contributor        |

## API

``` mermaid
flowchart LR
  file_dataset(Dataset)
  file_solution(Test data)
  file_embedding(Embedding)
  file_score(Score)
  file_common_dataset(Common dataset)
  comp_control_method[/Control method/]
  comp_method[/Method/]
  comp_metric[/Metric/]
  comp_process_dataset[/Data processor/]
  file_dataset---comp_control_method
  file_solution---comp_control_method
  file_dataset---comp_method
  file_embedding---comp_metric
  file_solution---comp_metric
  file_common_dataset---comp_process_dataset
  comp_control_method-->file_embedding
  comp_method-->file_embedding
  comp_metric-->file_score
  comp_process_dataset-->file_dataset
  comp_process_dataset-->file_solution
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
[`src/dimensionality_reduction`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/dimensionality_reduction)

A dimensionality reduction dataset processor.

Arguments:

<div class="small">

| Name                | Type   | Description                                                    |
|:--------------------|:-------|:---------------------------------------------------------------|
| `--input`           | `file` | A dataset processed by the common dataset processing pipeline. |
| `--output_dataset`  | `file` | (*Output*) The dataset to pass to a method.                    |
| `--output_solution` | `file` | (*Output*) The data for evaluating a dimensionality reduction. |

</div>

## File format: Dataset

The dataset to pass to a method.

Example file:
`resources_test/dimensionality_reduction/pancreas/dataset.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     var: 'hvg_score'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'normalization_id'

</div>

Slot description:

<div class="small">

| Slot                      | Type      | Description                                                                          |
|:--------------------------|:----------|:-------------------------------------------------------------------------------------|
| `var["hvg_score"]`        | `double`  | High variability gene score (normalized dispersion). The greater, the more variable. |
| `layers["counts"]`        | `integer` | Raw counts.                                                                          |
| `layers["normalized"]`    | `double`  | Normalized expression values.                                                        |
| `uns["dataset_id"]`       | `string`  | A unique identifier for the dataset.                                                 |
| `uns["normalization_id"]` | `string`  | Which normalization was used.                                                        |

</div>

## File format: Test data

The data for evaluating a dimensionality reduction.

Example file:
`resources_test/dimensionality_reduction/pancreas/solution.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     var: 'hvg_score'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'normalization_id'

</div>

Slot description:

<div class="small">

| Slot                      | Type      | Description                                                                          |
|:--------------------------|:----------|:-------------------------------------------------------------------------------------|
| `var["hvg_score"]`        | `double`  | High variability gene score (normalized dispersion). The greater, the more variable. |
| `layers["counts"]`        | `integer` | Raw counts.                                                                          |
| `layers["normalized"]`    | `double`  | Normalized expression values.                                                        |
| `uns["dataset_id"]`       | `string`  | A unique identifier for the dataset.                                                 |
| `uns["normalization_id"]` | `string`  | Which normalization was used.                                                        |

</div>

## Component type: Control method

Path:
[`src/dimensionality_reduction/control_methods`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/dimensionality_reduction/control_methods)

Quality control methods for verifying the pipeline.

Arguments:

<div class="small">

| Name               | Type   | Description                                                   |
|:-------------------|:-------|:--------------------------------------------------------------|
| `--input`          | `file` | The dataset to pass to a method.                              |
| `--input_solution` | `file` | The data for evaluating a dimensionality reduction.           |
| `--output`         | `file` | (*Output*) A dataset with dimensionality reduction embedding. |

</div>

## Component type: Method

Path:
[`src/dimensionality_reduction/methods`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/dimensionality_reduction/methods)

A dimensionality reduction method.

Arguments:

<div class="small">

| Name       | Type   | Description                                                   |
|:-----------|:-------|:--------------------------------------------------------------|
| `--input`  | `file` | The dataset to pass to a method.                              |
| `--output` | `file` | (*Output*) A dataset with dimensionality reduction embedding. |

</div>

## Component type: Metric

Path:
[`src/dimensionality_reduction/metrics`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/dimensionality_reduction/metrics)

A dimensionality reduction metric.

Arguments:

<div class="small">

| Name                | Type   | Description                                         |
|:--------------------|:-------|:----------------------------------------------------|
| `--input_embedding` | `file` | A dataset with dimensionality reduction embedding.  |
| `--input_solution`  | `file` | The data for evaluating a dimensionality reduction. |
| `--output`          | `file` | (*Output*) Metric score file.                       |

</div>

## File format: Embedding

A dataset with dimensionality reduction embedding.

Example file:
`resources_test/dimensionality_reduction/pancreas/embedding.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     obsm: 'X_emb'
     uns: 'dataset_id', 'method_id', 'normalization_id'

</div>

Slot description:

<div class="small">

| Slot                      | Type     | Description                          |
|:--------------------------|:---------|:-------------------------------------|
| `obsm["X_emb"]`           | `double` | The dimensionally reduced embedding. |
| `uns["dataset_id"]`       | `string` | A unique identifier for the dataset. |
| `uns["method_id"]`        | `string` | A unique identifier for the method.  |
| `uns["normalization_id"]` | `string` | Which normalization was used.        |

</div>

## File format: Score

Metric score file

Example file:
`resources_test/dimensionality_reduction/pancreas/score.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     uns: 'dataset_id', 'normalization_id', 'method_id', 'metric_ids', 'metric_values'

</div>

Slot description:

<div class="small">

| Slot                      | Type     | Description                                                                                  |
|:--------------------------|:---------|:---------------------------------------------------------------------------------------------|
| `uns["dataset_id"]`       | `string` | A unique identifier for the dataset.                                                         |
| `uns["normalization_id"]` | `string` | Which normalization was used.                                                                |
| `uns["method_id"]`        | `string` | A unique identifier for the method.                                                          |
| `uns["metric_ids"]`       | `string` | One or more unique metric identifiers.                                                       |
| `uns["metric_values"]`    | `double` | The metric values obtained for the given prediction. Must be of same length as ‘metric_ids’. |

</div>
