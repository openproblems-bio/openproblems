
# Label Projection

## The task

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

## The metrics

Metrics for label projection aim to characterize how well each classifer
correctly assigns cell type labels to cells in the test set.

- **Accuracy**: Average number of correctly applied labels.
- **F1 score**: The [F1
  score](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.f1_score.html)
  is a weighted average of the precision and recall over all class
  labels, where an F1 score reaches its best value at 1 and worst score
  at 0, where each class contributes to the score relative to its
  frequency in the dataset.
- **Macro F1 score**: The macro F1 score is an unweighted F1 score,
  where each class contributes equally, regardless of its frequency.

## API

``` mermaid
%%| column: screen-inset-shaded
flowchart LR
  dataset_censoring__output_train(Training data)
  dataset_censoring__output_test(Test data)
  dataset_censoring__output_solution(Solution)
  dataset_preprocessing__input(Raw dataset)
  dataset_preprocessing__output(Pre-processed dataset)
  method__output(Prediction)
  metric__output(Scores)
  dataset_censoring[/Dataset censoring/]
  dataset_preprocessing[/Dataset preprocessing/]
  method[/Method/]
  metric[/Metric/]
  dataset_preprocessing__output---dataset_censoring
  dataset_preprocessing__input---dataset_preprocessing
  dataset_censoring__output_train---method
  dataset_censoring__output_test---method
  dataset_censoring__output_solution---metric
  method__output---metric
  dataset_censoring-->dataset_censoring__output_train
  dataset_censoring-->dataset_censoring__output_test
  dataset_censoring-->dataset_censoring__output_solution
  dataset_preprocessing-->dataset_preprocessing__output
  method-->method__output
  metric-->metric__output
```

### Training data

The training data

Used in:

- dataset_censoring: output_train (as output)
- method: input_train (as input)

Slots:

| struct | name           | type    | description                                                         |
|:-------|:---------------|:--------|:--------------------------------------------------------------------|
| layers | counts         | integer | Raw counts                                                          |
| layers | lognorm        | double  | Log-transformed normalised counts                                   |
| obs    | labels         | double  | Ground truth cell type labels                                       |
| obs    | batch          | double  | Batch information                                                   |
| uns    | raw_dataset_id | string  | A unique identifier for the original dataset (before preprocessing) |
| uns    | dataset        | string  | A unique identifier for the dataset                                 |

### Test data

The censored test data

Used in:

- dataset_censoring: output_test (as output)
- method: input_test (as input)

Slots:

| struct | name           | type    | description                                                         |
|:-------|:---------------|:--------|:--------------------------------------------------------------------|
| layers | counts         | integer | Raw counts                                                          |
| layers | lognorm        | double  | Log-transformed normalised counts                                   |
| obs    | batch          | double  | Batch information                                                   |
| uns    | raw_dataset_id | string  | A unique identifier for the original dataset (before preprocessing) |
| uns    | dataset        | string  | A unique identifier for the dataset                                 |

### Solution

The solution for the test data

Used in:

- dataset_censoring: output_solution (as output)
- metric: input_solution (as input)

Slots:

| struct | name           | type    | description                                                         |
|:-------|:---------------|:--------|:--------------------------------------------------------------------|
| layers | counts         | integer | Raw counts                                                          |
| layers | lognorm        | double  | Log-transformed normalised counts                                   |
| obs    | labels         | double  | Ground truth cell type labels                                       |
| obs    | batch          | double  | Batch information                                                   |
| uns    | raw_dataset_id | string  | A unique identifier for the original dataset (before preprocessing) |
| uns    | dataset        | string  | A unique identifier for the dataset                                 |

### Raw dataset

An unprocessed dataset.

Used in:

- dataset_preprocessing: input (as input)

Slots:

| struct | name           | type    | description                                                         |
|:-------|:---------------|:--------|:--------------------------------------------------------------------|
| layers | counts         | integer | Raw counts                                                          |
| obs    | labels         | double  | Ground truth cell type labels                                       |
| obs    | batch          | double  | Batch information                                                   |
| uns    | raw_dataset_id | string  | A unique identifier for the original dataset (before preprocessing) |

### Pre-processed dataset

A preprocessed dataset

Used in:

- dataset_censoring: input (as input)
- dataset_preprocessing: output (as output)

Slots:

| struct | name           | type    | description                                                         |
|:-------|:---------------|:--------|:--------------------------------------------------------------------|
| layers | counts         | integer | Raw counts                                                          |
| layers | lognorm        | double  | Log-transformed normalised counts                                   |
| obs    | labels         | double  | Ground truth cell type labels                                       |
| obs    | batch          | double  | Batch information                                                   |
| uns    | raw_dataset_id | string  | A unique identifier for the original dataset (before preprocessing) |
| uns    | dataset        | string  | A unique identifier for the dataset                                 |

### Prediction

The prediction file

Used in:

- method: output (as output)
- metric: input_prediction (as input)

Slots:

| struct | name           | type   | description                                                         |
|:-------|:---------------|:-------|:--------------------------------------------------------------------|
| obs    | labels_pred    | double | Predicted labels for the test cells.                                |
| uns    | dataset_id     | string | A unique identifier for the dataset                                 |
| uns    | raw_dataset_id | string | A unique identifier for the original dataset (before preprocessing) |
| uns    | method_id      | string | A unique identifier for the method                                  |

### Scores

Metric score file

Used in:

- metric: output (as output)

Slots:

| struct | name           | type   | description                                                                                  |
|:-------|:---------------|:-------|:---------------------------------------------------------------------------------------------|
| uns    | dataset_id     | string | A unique identifier for the dataset                                                          |
| uns    | raw_dataset_id | string | A unique identifier for the original dataset (before preprocessing)                          |
| uns    | method_id      | string | A unique identifier for the method                                                           |
| uns    | metric_ids     | string | One or more unique metric identifiers                                                        |
| uns    | metric_values  | double | The metric values obtained for the given prediction. Must be of same length as ‘metric_ids’. |

<!--
Datasets should contain the following attributes:

* `adata.obs["labels"]` with ground truth celltype labels,
* `adata.obs["batch"]` with information of batches in the data, and
* `adata.obs["is_train"]` with a train vs. test split

It should be noted that datasets may only contain a single batch, or not contain discriminative batch information.

Methods should assign output celltype labels to `adata.obs['labels_pred']` using only the labels from the training data.

Note that the true labels are contained in `adata['labels']`.

Metrics can compare `adata['labels']` to `adata.obs['labels_pred']` using only the labels from the test data.
-->
