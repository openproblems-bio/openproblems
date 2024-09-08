# Label projection


Automated cell type annotation from rich, labeled reference data

Path:
[`src/tasks/label_projection`](https://github.com/openproblems-bio/openproblems/tree/main/src/tasks/label_projection)

## Motivation

A major challenge for integrating single cell datasets is creating
matching cell type annotations for each cell. One of the most common
strategies for annotating cell types is referred to as
[“cluster-then-annotate”](https://www.nature.com/articles/s41576-018-0088-9)
whereby cells are aggregated into clusters based on feature similarity
and then manually characterized based on differential gene expression or
previously identified marker genes. Recently, methods have emerged to
build on this strategy and annotate cells using [known marker
genes](https://www.nature.com/articles/s41592-019-0535-3). However,
these strategies pose a difficulty for integrating atlas-scale datasets
as the particular annotations may not match.

## Description

To ensure that the cell type labels in newly generated datasets match
existing reference datasets, some methods align cells to a previously
annotated [reference
dataset](https://academic.oup.com/bioinformatics/article/35/22/4688/54802990)
and then *project* labels from the reference to the new dataset.

Here, we compare methods for annotation based on a reference dataset.
The datasets consist of two or more samples of single cell profiles that
have been manually annotated with matching labels. These datasets are
then split into training and test batches, and the task of each method
is to train a cell type classifer on the training set and project those
labels onto the test set.

## Authors & contributors

| name              | roles              |
|:------------------|:-------------------|
| Nikolay Markov    | author, maintainer |
| Scott Gigante     | author             |
| Robrecht Cannoodt | author             |

## API

``` mermaid
flowchart LR
  file_common_dataset("Common Dataset")
  comp_process_dataset[/"Data processor"/]
  file_train("Training data")
  file_test("Test data")
  file_solution("Solution")
  comp_control_method[/"Control method"/]
  comp_method[/"Method"/]
  comp_metric[/"Metric"/]
  file_prediction("Prediction")
  file_score("Score")
  file_common_dataset---comp_process_dataset
  comp_process_dataset-->file_train
  comp_process_dataset-->file_test
  comp_process_dataset-->file_solution
  file_train---comp_control_method
  file_train---comp_method
  file_test---comp_control_method
  file_test---comp_method
  file_solution---comp_control_method
  file_solution---comp_metric
  comp_control_method-->file_prediction
  comp_method-->file_prediction
  comp_metric-->file_score
  file_prediction---comp_metric
```

## File format: Common Dataset

A subset of the common dataset.

Example file: `resources_test/common/pancreas/dataset.h5ad`

Format:

<div class="small">

    AnnData object
     obs: 'cell_type', 'batch'
     var: 'hvg', 'hvg_score'
     obsm: 'X_pca'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'dataset_name', 'dataset_url', 'dataset_reference', 'dataset_summary', 'dataset_description', 'dataset_organism', 'normalization_id'

</div>

Slot description:

<div class="small">

| Slot                         | Type      | Description                                                                    |
|:-----------------------------|:----------|:-------------------------------------------------------------------------------|
| `obs["cell_type"]`           | `string`  | Cell type information.                                                         |
| `obs["batch"]`               | `string`  | Batch information.                                                             |
| `var["hvg"]`                 | `boolean` | Whether or not the feature is considered to be a ‘highly variable gene’.       |
| `var["hvg_score"]`           | `double`  | A ranking of the features by hvg.                                              |
| `obsm["X_pca"]`              | `double`  | The resulting PCA embedding.                                                   |
| `layers["counts"]`           | `integer` | Raw counts.                                                                    |
| `layers["normalized"]`       | `double`  | Normalized expression values.                                                  |
| `uns["dataset_id"]`          | `string`  | A unique identifier for the dataset.                                           |
| `uns["dataset_name"]`        | `string`  | Nicely formatted name.                                                         |
| `uns["dataset_url"]`         | `string`  | (*Optional*) Link to the original source of the dataset.                       |
| `uns["dataset_reference"]`   | `string`  | (*Optional*) Bibtex reference of the paper in which the dataset was published. |
| `uns["dataset_summary"]`     | `string`  | Short description of the dataset.                                              |
| `uns["dataset_description"]` | `string`  | Long description of the dataset.                                               |
| `uns["dataset_organism"]`    | `string`  | (*Optional*) The organism of the sample in the dataset.                        |
| `uns["normalization_id"]`    | `string`  | Which normalization was used.                                                  |

</div>

## Component type: Data processor

Path:
[`src/label_projection`](https://github.com/openproblems-bio/openproblems/tree/main/src/label_projection)

A label projection dataset processor.

Arguments:

<div class="small">

| Name                | Type   | Description                                |
|:--------------------|:-------|:-------------------------------------------|
| `--input`           | `file` | A subset of the common dataset.            |
| `--output_train`    | `file` | (*Output*) The training data.              |
| `--output_test`     | `file` | (*Output*) The test data (without labels). |
| `--output_solution` | `file` | (*Output*) The solution for the test data. |

</div>

## File format: Training data

The training data

Example file: `resources_test/label_projection/pancreas/train.h5ad`

Format:

<div class="small">

    AnnData object
     obs: 'label', 'batch'
     var: 'hvg', 'hvg_score'
     obsm: 'X_pca'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'normalization_id'

</div>

Slot description:

<div class="small">

| Slot                      | Type      | Description                                                              |
|:--------------------------|:----------|:-------------------------------------------------------------------------|
| `obs["label"]`            | `string`  | Ground truth cell type labels.                                           |
| `obs["batch"]`            | `string`  | Batch information.                                                       |
| `var["hvg"]`              | `boolean` | Whether or not the feature is considered to be a ‘highly variable gene’. |
| `var["hvg_score"]`        | `double`  | A ranking of the features by hvg.                                        |
| `obsm["X_pca"]`           | `double`  | The resulting PCA embedding.                                             |
| `layers["counts"]`        | `integer` | Raw counts.                                                              |
| `layers["normalized"]`    | `double`  | Normalized counts.                                                       |
| `uns["dataset_id"]`       | `string`  | A unique identifier for the dataset.                                     |
| `uns["normalization_id"]` | `string`  | Which normalization was used.                                            |

</div>

## File format: Test data

The test data (without labels)

Example file: `resources_test/label_projection/pancreas/test.h5ad`

Format:

<div class="small">

    AnnData object
     obs: 'batch'
     var: 'hvg', 'hvg_score'
     obsm: 'X_pca'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'normalization_id'

</div>

Slot description:

<div class="small">

| Slot                      | Type      | Description                                                              |
|:--------------------------|:----------|:-------------------------------------------------------------------------|
| `obs["batch"]`            | `string`  | Batch information.                                                       |
| `var["hvg"]`              | `boolean` | Whether or not the feature is considered to be a ‘highly variable gene’. |
| `var["hvg_score"]`        | `double`  | A ranking of the features by hvg.                                        |
| `obsm["X_pca"]`           | `double`  | The resulting PCA embedding.                                             |
| `layers["counts"]`        | `integer` | Raw counts.                                                              |
| `layers["normalized"]`    | `double`  | Normalized counts.                                                       |
| `uns["dataset_id"]`       | `string`  | A unique identifier for the dataset.                                     |
| `uns["normalization_id"]` | `string`  | Which normalization was used.                                            |

</div>

## File format: Solution

The solution for the test data

Example file: `resources_test/label_projection/pancreas/solution.h5ad`

Format:

<div class="small">

    AnnData object
     obs: 'label', 'batch'
     var: 'hvg', 'hvg_score'
     obsm: 'X_pca'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'dataset_name', 'dataset_url', 'dataset_reference', 'dataset_summary', 'dataset_description', 'dataset_organism', 'normalization_id'

</div>

Slot description:

<div class="small">

| Slot                         | Type      | Description                                                                    |
|:-----------------------------|:----------|:-------------------------------------------------------------------------------|
| `obs["label"]`               | `string`  | Ground truth cell type labels.                                                 |
| `obs["batch"]`               | `string`  | Batch information.                                                             |
| `var["hvg"]`                 | `boolean` | Whether or not the feature is considered to be a ‘highly variable gene’.       |
| `var["hvg_score"]`           | `double`  | A ranking of the features by hvg.                                              |
| `obsm["X_pca"]`              | `double`  | The resulting PCA embedding.                                                   |
| `layers["counts"]`           | `integer` | Raw counts.                                                                    |
| `layers["normalized"]`       | `double`  | Normalized counts.                                                             |
| `uns["dataset_id"]`          | `string`  | A unique identifier for the dataset.                                           |
| `uns["dataset_name"]`        | `string`  | Nicely formatted name.                                                         |
| `uns["dataset_url"]`         | `string`  | (*Optional*) Link to the original source of the dataset.                       |
| `uns["dataset_reference"]`   | `string`  | (*Optional*) Bibtex reference of the paper in which the dataset was published. |
| `uns["dataset_summary"]`     | `string`  | Short description of the dataset.                                              |
| `uns["dataset_description"]` | `string`  | Long description of the dataset.                                               |
| `uns["dataset_organism"]`    | `string`  | (*Optional*) The organism of the sample in the dataset.                        |
| `uns["normalization_id"]`    | `string`  | Which normalization was used.                                                  |

</div>

## Component type: Control method

Path:
[`src/label_projection/control_methods`](https://github.com/openproblems-bio/openproblems/tree/main/src/label_projection/control_methods)

Quality control methods for verifying the pipeline.

Arguments:

<div class="small">

| Name               | Type   | Description                     |
|:-------------------|:-------|:--------------------------------|
| `--input_train`    | `file` | The training data.              |
| `--input_test`     | `file` | The test data (without labels). |
| `--input_solution` | `file` | The solution for the test data. |
| `--output`         | `file` | (*Output*) The prediction file. |

</div>

## Component type: Method

Path:
[`src/label_projection/methods`](https://github.com/openproblems-bio/openproblems/tree/main/src/label_projection/methods)

A label projection method.

Arguments:

<div class="small">

| Name            | Type   | Description                     |
|:----------------|:-------|:--------------------------------|
| `--input_train` | `file` | The training data.              |
| `--input_test`  | `file` | The test data (without labels). |
| `--output`      | `file` | (*Output*) The prediction file. |

</div>

## Component type: Metric

Path:
[`src/label_projection/metrics`](https://github.com/openproblems-bio/openproblems/tree/main/src/label_projection/metrics)

A label projection metric.

Arguments:

<div class="small">

| Name                 | Type   | Description                     |
|:---------------------|:-------|:--------------------------------|
| `--input_solution`   | `file` | The solution for the test data. |
| `--input_prediction` | `file` | The prediction file.            |
| `--output`           | `file` | (*Output*) Metric score file.   |

</div>

## File format: Prediction

The prediction file

Example file: `resources_test/label_projection/pancreas/prediction.h5ad`

Format:

<div class="small">

    AnnData object
     obs: 'label_pred'
     uns: 'dataset_id', 'normalization_id', 'method_id'

</div>

Slot description:

<div class="small">

| Slot                      | Type     | Description                          |
|:--------------------------|:---------|:-------------------------------------|
| `obs["label_pred"]`       | `string` | Predicted labels for the test cells. |
| `uns["dataset_id"]`       | `string` | A unique identifier for the dataset. |
| `uns["normalization_id"]` | `string` | Which normalization was used.        |
| `uns["method_id"]`        | `string` | A unique identifier for the method.  |

</div>

## File format: Score

Metric score file

Example file: `resources_test/label_projection/pancreas/score.h5ad`

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

