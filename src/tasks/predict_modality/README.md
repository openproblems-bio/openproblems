# Predict Modality

Predicting the profiles of one modality (e.g. protein abundance) from
another (e.g. mRNA expression).

Path:
[`src/tasks/predict_modality`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/tasks/predict_modality)

## Motivation

Experimental techniques to measure multiple modalities within the same
single cell are increasingly becoming available. The demand for these
measurements is driven by the promise to provide a deeper insight into
the state of a cell. Yet, the modalities are also intrinsically linked.
We know that DNA must be accessible (ATAC data) to produce mRNA
(expression data), and mRNA in turn is used as a template to produce
protein (protein abundance). These processes are regulated often by the
same molecules that they produce: for example, a protein may bind DNA to
prevent the production of more mRNA. Understanding these regulatory
processes would be transformative for synthetic biology and drug target
discovery. Any method that can predict a modality from another must have
accounted for these regulatory processes, but the demand for multi-modal
data shows that this is not trivial.

## Description

In this task, the goal is to take one modality and predict the other
modality for all features in each cell. This task requires translating
information between multiple layers of gene regulation. In some ways,
this is similar to the task of machine translation. In machine
translation, the same sentiment is expressed in multiple languages and
the goal is to train a model to represent the same meaning in a
different language. In this context, the same cellular state is measured
in two different feature sets and the goal of this task is to translate
the information about cellular state from one modality to the other.

## Authors & contributors

| name               | roles              |
|:-------------------|:-------------------|
| Robrecht Cannoodt  | author, maintainer |
| Kai Waldrant       | contributor        |
| Louise Deconinck   | author             |
| Alex Tong          | author             |
| Bastian Rieck      | author             |
| Daniel Burkhardt   | author             |
| Alejandro Granados | author             |

## API

``` mermaid
flowchart LR
  file_dataset_rna("Raw dataset RNA")
  comp_process_dataset[/"Data processor"/]
  file_train_mod1("Train mod1")
  file_train_mod2("Train mod2")
  file_test_mod1("Test mod1")
  file_test_mod2("Test mod2")
  comp_control_method[/"Control method"/]
  comp_method[/"Method"/]
  comp_metric[/"Metric"/]
  file_prediction("Prediction")
  file_score("Score")
  file_dataset_other_mod("Raw dataset mod2")
  file_dataset_rna---comp_process_dataset
  comp_process_dataset-->file_train_mod1
  comp_process_dataset-->file_train_mod2
  comp_process_dataset-->file_test_mod1
  comp_process_dataset-->file_test_mod2
  file_train_mod1---comp_control_method
  file_train_mod1---comp_method
  file_train_mod2---comp_control_method
  file_train_mod2---comp_method
  file_test_mod1---comp_control_method
  file_test_mod1---comp_method
  file_test_mod2---comp_control_method
  file_test_mod2---comp_metric
  comp_control_method-->file_prediction
  comp_method-->file_prediction
  comp_metric-->file_score
  file_prediction---comp_metric
  file_dataset_other_mod---comp_process_dataset
```

## File format: Raw dataset RNA

The RNA modality of the raw dataset.

Example file:
`resources_test/common/bmmc_cite_starter/openproblems_bmmc_cite_starter.output_rna.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     obs: 'batch', 'size_factors'
     var: 'gene_ids'
     obsm: 'gene_activity'
     layers: 'counts'
     uns: 'dataset_id', 'gene_activity_var_names'

</div>

Slot description:

<div class="small">

| Slot                             | Type      | Description                                                        |
|:---------------------------------|:----------|:-------------------------------------------------------------------|
| `obs["batch"]`                   | `string`  | Batch information.                                                 |
| `obs["size_factors"]`            | `double`  | (*Optional*) The size factors of the cells prior to normalization. |
| `var["gene_ids"]`                | `string`  | (*Optional*) The gene identifiers (if available).                  |
| `obsm["gene_activity"]`          | `double`  | (*Optional*) ATAC gene activity.                                   |
| `layers["counts"]`               | `integer` | Raw counts.                                                        |
| `uns["dataset_id"]`              | `string`  | A unique identifier for the dataset.                               |
| `uns["gene_activity_var_names"]` | `string`  | (*Optional*) Names of the gene activity matrix.                    |

</div>

## Component type: Data processor

Path:
[`src/predict_modality`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/predict_modality)

A predict modality dataset processor.

Arguments:

<div class="small">

| Name                  | Type      | Description                                                                |
|:----------------------|:----------|:---------------------------------------------------------------------------|
| `--input_rna`         | `file`    | The RNA modality of the raw dataset.                                       |
| `--input_other_mod`   | `file`    | The second modality of the raw dataset. Must be an ADT or an ATAC dataset. |
| `--output_train_mod1` | `file`    | (*Output*) The mod1 expression values of the train cells.                  |
| `--output_train_mod2` | `file`    | (*Output*) The mod2 expression values of the train cells.                  |
| `--output_test_mod1`  | `file`    | (*Output*) The mod1 expression values of the test cells.                   |
| `--output_test_mod2`  | `file`    | (*Output*) The mod2 expression values of the test cells.                   |
| `--seed`              | `integer` | (*Optional*) The seed for determining the train/test split. Default: `1`.  |

</div>

## File format: Train mod1

The mod1 expression values of the train cells.

Example file:
`resources_test/predict_modality/bmmc_cite_starter/train_mod1.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     obs: 'batch', 'size_factors'
     var: 'gene_ids'
     obsm: 'gene_activity'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'dataset_organism', 'gene_activity_var_names'

</div>

Slot description:

<div class="small">

| Slot                             | Type      | Description                                                        |
|:---------------------------------|:----------|:-------------------------------------------------------------------|
| `obs["batch"]`                   | `string`  | Batch information.                                                 |
| `obs["size_factors"]`            | `double`  | (*Optional*) The size factors of the cells prior to normalization. |
| `var["gene_ids"]`                | `string`  | (*Optional*) The gene identifiers (if available).                  |
| `obsm["gene_activity"]`          | `double`  | (*Optional*) ATAC gene activity.                                   |
| `layers["counts"]`               | `integer` | Raw counts.                                                        |
| `layers["normalized"]`           | `double`  | Normalized expression values.                                      |
| `uns["dataset_id"]`              | `string`  | A unique identifier for the dataset.                               |
| `uns["dataset_organism"]`        | `string`  | (*Optional*) The organism of the sample in the dataset.            |
| `uns["gene_activity_var_names"]` | `string`  | (*Optional*) Names of the gene activity matrix.                    |

</div>

## File format: Train mod2

The mod2 expression values of the train cells.

Example file:
`resources_test/predict_modality/bmmc_cite_starter/train_mod2.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     obs: 'batch', 'size_factors'
     var: 'gene_ids'
     obsm: 'gene_activity'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'dataset_organism', 'gene_activity_var_names'

</div>

Slot description:

<div class="small">

| Slot                             | Type      | Description                                                        |
|:---------------------------------|:----------|:-------------------------------------------------------------------|
| `obs["batch"]`                   | `string`  | Batch information.                                                 |
| `obs["size_factors"]`            | `double`  | (*Optional*) The size factors of the cells prior to normalization. |
| `var["gene_ids"]`                | `string`  | (*Optional*) The gene identifiers (if available).                  |
| `obsm["gene_activity"]`          | `double`  | (*Optional*) ATAC gene activity.                                   |
| `layers["counts"]`               | `integer` | Raw counts.                                                        |
| `layers["normalized"]`           | `double`  | Normalized expression values.                                      |
| `uns["dataset_id"]`              | `string`  | A unique identifier for the dataset.                               |
| `uns["dataset_organism"]`        | `string`  | (*Optional*) The organism of the sample in the dataset.            |
| `uns["gene_activity_var_names"]` | `string`  | (*Optional*) Names of the gene activity matrix.                    |

</div>

## File format: Test mod1

The mod1 expression values of the test cells.

Example file:
`resources_test/predict_modality/bmmc_cite_starter/test_mod1.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     obs: 'batch', 'size_factors'
     var: 'gene_ids'
     obsm: 'gene_activity'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'dataset_organism', 'gene_activity_var_names'

</div>

Slot description:

<div class="small">

| Slot                             | Type      | Description                                                        |
|:---------------------------------|:----------|:-------------------------------------------------------------------|
| `obs["batch"]`                   | `string`  | Batch information.                                                 |
| `obs["size_factors"]`            | `double`  | (*Optional*) The size factors of the cells prior to normalization. |
| `var["gene_ids"]`                | `string`  | (*Optional*) The gene identifiers (if available).                  |
| `obsm["gene_activity"]`          | `double`  | (*Optional*) ATAC gene activity.                                   |
| `layers["counts"]`               | `integer` | Raw counts.                                                        |
| `layers["normalized"]`           | `double`  | Normalized expression values.                                      |
| `uns["dataset_id"]`              | `string`  | A unique identifier for the dataset.                               |
| `uns["dataset_organism"]`        | `string`  | (*Optional*) The organism of the sample in the dataset.            |
| `uns["gene_activity_var_names"]` | `string`  | (*Optional*) Names of the gene activity matrix.                    |

</div>

## File format: Test mod2

The mod2 expression values of the test cells.

Example file:
`resources_test/predict_modality/bmmc_cite_starter/test_mod2.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     obs: 'batch', 'size_factors'
     var: 'gene_ids'
     obsm: 'gene_activity'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'dataset_organism', 'gene_activity_var_names'

</div>

Slot description:

<div class="small">

| Slot                             | Type      | Description                                                        |
|:---------------------------------|:----------|:-------------------------------------------------------------------|
| `obs["batch"]`                   | `string`  | Batch information.                                                 |
| `obs["size_factors"]`            | `double`  | (*Optional*) The size factors of the cells prior to normalization. |
| `var["gene_ids"]`                | `string`  | (*Optional*) The gene identifiers (if available).                  |
| `obsm["gene_activity"]`          | `double`  | (*Optional*) ATAC gene activity.                                   |
| `layers["counts"]`               | `integer` | Raw counts.                                                        |
| `layers["normalized"]`           | `double`  | Normalized expression values.                                      |
| `uns["dataset_id"]`              | `string`  | A unique identifier for the dataset.                               |
| `uns["dataset_organism"]`        | `string`  | (*Optional*) The organism of the sample in the dataset.            |
| `uns["gene_activity_var_names"]` | `string`  | (*Optional*) Names of the gene activity matrix.                    |

</div>

## Component type: Control method

Path:
[`src/predict_modality/control_methods`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/predict_modality/control_methods)

Quality control methods for verifying the pipeline.

Arguments:

<div class="small">

| Name                 | Type   | Description                                                              |
|:---------------------|:-------|:-------------------------------------------------------------------------|
| `--input_train_mod1` | `file` | The mod1 expression values of the train cells.                           |
| `--input_train_mod2` | `file` | The mod2 expression values of the train cells.                           |
| `--input_test_mod1`  | `file` | The mod1 expression values of the test cells.                            |
| `--input_test_mod2`  | `file` | The mod2 expression values of the test cells.                            |
| `--output`           | `file` | (*Output*) A prediction of the mod2 expression values of the test cells. |

</div>

## Component type: Method

Path:
[`src/predict_modality/methods`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/predict_modality/methods)

A regression method.

Arguments:

<div class="small">

| Name                 | Type   | Description                                                              |
|:---------------------|:-------|:-------------------------------------------------------------------------|
| `--input_train_mod1` | `file` | The mod1 expression values of the train cells.                           |
| `--input_train_mod2` | `file` | The mod2 expression values of the train cells.                           |
| `--input_test_mod1`  | `file` | The mod1 expression values of the test cells.                            |
| `--output`           | `file` | (*Output*) A prediction of the mod2 expression values of the test cells. |

</div>

## Component type: Metric

Path:
[`src/predict_modality/metrics`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/predict_modality/metrics)

A predict modality metric.

Arguments:

<div class="small">

| Name                 | Type   | Description                                                   |
|:---------------------|:-------|:--------------------------------------------------------------|
| `--input_prediction` | `file` | A prediction of the mod2 expression values of the test cells. |
| `--input_test_mod2`  | `file` | The mod2 expression values of the test cells.                 |
| `--output`           | `file` | (*Output*) Metric score file.                                 |

</div>

## File format: Prediction

A prediction of the mod2 expression values of the test cells

Example file:
`resources_test/predict_modality/bmmc_cite_starter/prediction.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     layers: 'normalized'
     uns: 'dataset_id', 'method_id'

</div>

Slot description:

<div class="small">

| Slot                   | Type     | Description                             |
|:-----------------------|:---------|:----------------------------------------|
| `layers["normalized"]` | `double` | Predicted normalized expression values. |
| `uns["dataset_id"]`    | `string` | A unique identifier for the dataset.    |
| `uns["method_id"]`     | `string` | A unique identifier for the method.     |

</div>

## File format: Score

Metric score file

Example file:
`resources_test/predict_modality/bmmc_cite_starter/score.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     uns: 'dataset_id', 'method_id', 'metric_ids', 'metric_values'

</div>

Slot description:

<div class="small">

| Slot                   | Type     | Description                                                                                  |
|:-----------------------|:---------|:---------------------------------------------------------------------------------------------|
| `uns["dataset_id"]`    | `string` | A unique identifier for the dataset.                                                         |
| `uns["method_id"]`     | `string` | A unique identifier for the method.                                                          |
| `uns["metric_ids"]`    | `string` | One or more unique metric identifiers.                                                       |
| `uns["metric_values"]` | `double` | The metric values obtained for the given prediction. Must be of same length as ‘metric_ids’. |

</div>

## File format: Raw dataset mod2

The second modality of the raw dataset. Must be an ADT or an ATAC
dataset

Example file:
`resources_test/common/bmmc_cite_starter/openproblems_bmmc_cite_starter.output_adt.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     obs: 'batch', 'size_factors'
     var: 'gene_ids'
     obsm: 'gene_activity'
     layers: 'counts'
     uns: 'dataset_id', 'gene_activity_var_names'

</div>

Slot description:

<div class="small">

| Slot                             | Type      | Description                                                        |
|:---------------------------------|:----------|:-------------------------------------------------------------------|
| `obs["batch"]`                   | `string`  | Batch information.                                                 |
| `obs["size_factors"]`            | `double`  | (*Optional*) The size factors of the cells prior to normalization. |
| `var["gene_ids"]`                | `string`  | (*Optional*) The gene identifiers (if available).                  |
| `obsm["gene_activity"]`          | `double`  | (*Optional*) ATAC gene activity.                                   |
| `layers["counts"]`               | `integer` | Raw counts.                                                        |
| `uns["dataset_id"]`              | `string`  | A unique identifier for the dataset.                               |
| `uns["gene_activity_var_names"]` | `string`  | (*Optional*) Names of the gene activity matrix.                    |

</div>
