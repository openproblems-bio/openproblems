# Match Modalities

Match cells across datasets of the same set of samples on different
technologies / modalities.

Path:
[`src/tasks/match_modalities`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/tasks/match_modalities)

## Motivation

Cellular function is regulated by the complex interplay of different
types of biological molecules (DNA, RNA, proteins, etc.), which
determine the state of a cell. Several recently described technologies
allow for simultaneous measurement of different aspects of cellular
state. For example, sci-CAR \[@cao2018joint\] jointly profiles RNA
expression and chromatin accessibility on the same cell and CITE-seq
\[@stoeckius2017simultaneous\] measures surface protein abundance and
RNA expression from each cell. These technologies enable us to better
understand cellular function, however datasets are still rare and there
are tradeoffs that these measurements make for to profile multiple
modalities.

Joint methods can be more expensive or lower throughput or more noisy
than measuring a single modality at a time. Therefore it is useful to
develop methods that are capable of integrating measurements of the same
biological system but obtained using different technologies on different
cells.

## Description

In this task, the goal is to learn a latent space where cells profiled
by different technologies in different modalities are matched if they
have the same state. We use jointly profiled data as ground truth so
that we can evaluate when the observations from the same cell acquired
using different modalities are similar. A perfect result has each of the
paired observations sharing the same coordinates in the latent space. A
method that can achieve this would be able to match datasets across
modalities to enable multimodal cellular analysis from separately
measured profiles.

## Authors & contributors

| name              | roles              |
|:------------------|:-------------------|
| Scott Gigante     | author, maintainer |
| Alex Tong         | author             |
| Robrecht Cannoodt | author             |
| Kai Waldrant      | contributor        |

## API

``` mermaid
flowchart LR
  file_common_dataset_mod1("Common dataset mod1")
  comp_process_dataset[/"Data processor"/]
  file_dataset_mod1("Modality 1")
  file_dataset_mod2("Modality 2")
  file_solution_mod1("Solution mod1")
  file_solution_mod2("Solution mod1")
  comp_control_method[/"Control method"/]
  comp_method[/"Method"/]
  comp_metric[/"Metric"/]
  file_integrated_mod1("Integrated mod1")
  file_integrated_mod2("Integrated mod2")
  file_score("Score")
  file_common_dataset_mod2("Common dataset mod2")
  file_common_dataset_mod1---comp_process_dataset
  comp_process_dataset-->file_dataset_mod1
  comp_process_dataset-->file_dataset_mod2
  comp_process_dataset-->file_solution_mod1
  comp_process_dataset-->file_solution_mod2
  file_dataset_mod1---comp_control_method
  file_dataset_mod1---comp_method
  file_dataset_mod2---comp_control_method
  file_dataset_mod2---comp_method
  file_solution_mod1---comp_control_method
  file_solution_mod1---comp_metric
  file_solution_mod2---comp_control_method
  file_solution_mod2---comp_metric
  comp_control_method-->file_integrated_mod1
  comp_control_method-->file_integrated_mod2
  comp_method-->file_integrated_mod1
  comp_method-->file_integrated_mod2
  comp_metric-->file_score
  file_integrated_mod1---comp_metric
  file_integrated_mod2---comp_metric
  file_common_dataset_mod2---comp_process_dataset
```

## File format: Common dataset mod1

The first modality (RNA) of a dataset processed by the common multimodal
dataset processing pipeline.

Example file:
`resources_test/common/scicar_cell_lines/dataset_mod1.h5ad`

Description:

This dataset contains both raw counts and normalized data matrices, as
well as a PCA embedding, HVG selection and a kNN graph.

Format:

<div class="small">

    AnnData object
     obsm: 'X_svd'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'dataset_name', 'dataset_url', 'dataset_reference', 'dataset_summary', 'dataset_description', 'dataset_organism', 'normalization_id'

</div>

Slot description:

<div class="small">

| Slot                         | Type      | Description                                                                    |
|:-----------------------------|:----------|:-------------------------------------------------------------------------------|
| `obsm["X_svd"]`              | `double`  | The resulting SVD PCA embedding.                                               |
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

## Component type: Data processor

Path:
[`src/match_modalities`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/match_modalities)

A match modalities dataset processor.

Arguments:

<div class="small">

| Name                     | Type   | Description                                                                                                    |
|:-------------------------|:-------|:---------------------------------------------------------------------------------------------------------------|
| `--input_mod1`           | `file` | The first modality (RNA) of a dataset processed by the common multimodal dataset processing pipeline.          |
| `--input_mod2`           | `file` | The second modality (ADT or ATAC) of a dataset processed by the common multimodal dataset processing pipeline. |
| `--output_mod1`          | `file` | (*Output*) The first modality of a multimodal dataset. The cells of this dataset are randomly permuted.        |
| `--output_mod2`          | `file` | (*Output*) The second modality of a multimodal dataset. The cells of this dataset are randomly permuted.       |
| `--output_solution_mod1` | `file` | (*Output*) The ground truth information for the first modality.                                                |
| `--output_solution_mod2` | `file` | (*Output*) The ground truth information for the second modality.                                               |

</div>

## File format: Modality 1

The first modality of a multimodal dataset. The cells of this dataset
are randomly permuted.

Example file:
`resources_test/match_modalities/scicar_cell_lines/dataset_mod1.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     obsm: 'X_svd'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'normalization_id'

</div>

Slot description:

<div class="small">

| Slot                      | Type      | Description                          |
|:--------------------------|:----------|:-------------------------------------|
| `obsm["X_svd"]`           | `double`  | The resulting SVD PCA embedding.     |
| `layers["counts"]`        | `integer` | Raw counts.                          |
| `layers["normalized"]`    | `double`  | Normalized counts.                   |
| `uns["dataset_id"]`       | `string`  | A unique identifier for the dataset. |
| `uns["normalization_id"]` | `string`  | Which normalization was used.        |

</div>

## File format: Modality 2

The second modality of a multimodal dataset. The cells of this dataset
are randomly permuted.

Example file:
`resources_test/match_modalities/scicar_cell_lines/dataset_mod2.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     obsm: 'X_svd'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'normalization_id'

</div>

Slot description:

<div class="small">

| Slot                      | Type      | Description                          |
|:--------------------------|:----------|:-------------------------------------|
| `obsm["X_svd"]`           | `double`  | The resulting SVD PCA embedding.     |
| `layers["counts"]`        | `integer` | Raw counts.                          |
| `layers["normalized"]`    | `double`  | Normalized counts.                   |
| `uns["dataset_id"]`       | `string`  | A unique identifier for the dataset. |
| `uns["normalization_id"]` | `string`  | Which normalization was used.        |

</div>

## File format: Solution mod1

The ground truth information for the first modality

Example file:
`resources_test/match_modalities/scicar_cell_lines/solution_mod1.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     obs: 'permutation_indices'
     obsm: 'X_svd'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'dataset_name', 'dataset_url', 'dataset_reference', 'dataset_summary', 'dataset_description', 'dataset_organism', 'normalization_id'

</div>

Slot description:

<div class="small">

| Slot                         | Type      | Description                                                                    |
|:-----------------------------|:----------|:-------------------------------------------------------------------------------|
| `obs["permutation_indices"]` | `integer` | Indices with which to revert the permutation of the cells.                     |
| `obsm["X_svd"]`              | `double`  | The resulting SVD PCA embedding.                                               |
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

## File format: Solution mod1

The ground truth information for the second modality

Example file:
`resources_test/match_modalities/scicar_cell_lines/solution_mod2.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     obs: 'permutation_indices'
     obsm: 'X_svd'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'dataset_name', 'dataset_url', 'dataset_reference', 'dataset_summary', 'dataset_description', 'dataset_organism', 'normalization_id'

</div>

Slot description:

<div class="small">

| Slot                         | Type      | Description                                                                    |
|:-----------------------------|:----------|:-------------------------------------------------------------------------------|
| `obs["permutation_indices"]` | `integer` | Indices with which to revert the permutation of the cells.                     |
| `obsm["X_svd"]`              | `double`  | The resulting SVD PCA embedding.                                               |
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
[`src/match_modalities/control_methods`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/match_modalities/control_methods)

A multimodal data integration control method.

Arguments:

<div class="small">

| Name                    | Type   | Description                                                                                   |
|:------------------------|:-------|:----------------------------------------------------------------------------------------------|
| `--input_mod1`          | `file` | The first modality of a multimodal dataset. The cells of this dataset are randomly permuted.  |
| `--input_mod2`          | `file` | The second modality of a multimodal dataset. The cells of this dataset are randomly permuted. |
| `--input_solution_mod1` | `file` | The ground truth information for the first modality.                                          |
| `--input_solution_mod2` | `file` | The ground truth information for the second modality.                                         |
| `--output_mod1`         | `file` | (*Output*) The integrated embedding for the first modality.                                   |
| `--output_mod2`         | `file` | (*Output*) The integrated embedding for the second modality.                                  |

</div>

## Component type: Method

Path:
[`src/match_modalities/methods`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/match_modalities/methods)

A multimodal data integration method.

Arguments:

<div class="small">

| Name            | Type   | Description                                                                                   |
|:----------------|:-------|:----------------------------------------------------------------------------------------------|
| `--input_mod1`  | `file` | The first modality of a multimodal dataset. The cells of this dataset are randomly permuted.  |
| `--input_mod2`  | `file` | The second modality of a multimodal dataset. The cells of this dataset are randomly permuted. |
| `--output_mod1` | `file` | (*Output*) The integrated embedding for the first modality.                                   |
| `--output_mod2` | `file` | (*Output*) The integrated embedding for the second modality.                                  |

</div>

## Component type: Metric

Path:
[`src/match_modalities/metrics`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/match_modalities/metrics)

A multimodal data integration metric.

Arguments:

<div class="small">

| Name                      | Type   | Description                                           |
|:--------------------------|:-------|:------------------------------------------------------|
| `--input_integrated_mod1` | `file` | The integrated embedding for the first modality.      |
| `--input_integrated_mod2` | `file` | The integrated embedding for the second modality.     |
| `--input_solution_mod1`   | `file` | The ground truth information for the first modality.  |
| `--input_solution_mod2`   | `file` | The ground truth information for the second modality. |
| `--output`                | `file` | (*Output*) Metric score file.                         |

</div>

## File format: Integrated mod1

The integrated embedding for the first modality

Example file:
`resources_test/match_modalities/scicar_cell_lines/integrated_mod1.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     obsm: 'integrated'
     uns: 'dataset_id', 'normalization_id', 'method_id'

</div>

Slot description:

<div class="small">

| Slot                      | Type     | Description                          |
|:--------------------------|:---------|:-------------------------------------|
| `obsm["integrated"]`      | `double` | An integrated embedding.             |
| `uns["dataset_id"]`       | `string` | A unique identifier for the dataset. |
| `uns["normalization_id"]` | `string` | Which normalization was used.        |
| `uns["method_id"]`        | `string` | Which method was used.               |

</div>

## File format: Integrated mod2

The integrated embedding for the second modality

Example file:
`resources_test/match_modalities/scicar_cell_lines/integrated_mod2.h5ad`

Description:

NA

Format:

<div class="small">

    AnnData object
     obsm: 'integrated'
     uns: 'dataset_id', 'normalization_id', 'method_id'

</div>

Slot description:

<div class="small">

| Slot                      | Type     | Description                          |
|:--------------------------|:---------|:-------------------------------------|
| `obsm["integrated"]`      | `double` | An integrated embedding.             |
| `uns["dataset_id"]`       | `string` | A unique identifier for the dataset. |
| `uns["normalization_id"]` | `string` | Which normalization was used.        |
| `uns["method_id"]`        | `string` | Which method was used.               |

</div>

## File format: Score

Metric score file

Example file:
`resources_test/match_modalities/scicar_cell_lines/score.h5ad`

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

## File format: Common dataset mod2

The second modality (ADT or ATAC) of a dataset processed by the common
multimodal dataset processing pipeline.

Example file:
`resources_test/common/scicar_cell_lines/dataset_mod2.h5ad`

Description:

This dataset contains both raw counts and normalized data matrices, as
well as a PCA embedding, HVG selection and a kNN graph.

Format:

<div class="small">

    AnnData object
     obsm: 'X_svd'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'dataset_name', 'dataset_url', 'dataset_reference', 'dataset_summary', 'dataset_description', 'dataset_organism', 'normalization_id'

</div>

Slot description:

<div class="small">

| Slot                         | Type      | Description                                                                    |
|:-----------------------------|:----------|:-------------------------------------------------------------------------------|
| `obsm["X_svd"]`              | `double`  | The resulting SVD PCA embedding.                                               |
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
