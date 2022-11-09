
- <a href="#label-projection" id="toc-label-projection">Label
  Projection</a>
  - <a href="#the-task" id="toc-the-task">The task</a>
  - <a href="#methods" id="toc-methods">Methods</a>
  - <a href="#metrics" id="toc-metrics">Metrics</a>
  - <a href="#pipeline-topology" id="toc-pipeline-topology">Pipeline
    topology</a>
  - <a href="#file-format-api" id="toc-file-format-api">File format API</a>
    - <a href="#file-dataset" id="toc-file-dataset"><code>dataset.h5ad</code>:
      Preprocessed dataset</a>
    - <a href="#file-prediction"
      id="toc-file-prediction"><code>prediction.h5ad</code>: Prediction</a>
    - <a href="#file-score" id="toc-file-score"><code>score.h5ad</code>:
      Score</a>
    - <a href="#file-solution"
      id="toc-file-solution"><code>solution.h5ad</code>: Solution</a>
    - <a href="#file-test" id="toc-file-test"><code>test.h5ad</code>: Test
      data</a>
    - <a href="#file-train" id="toc-file-train"><code>train.h5ad</code>:
      Training data</a>
  - <a href="#component-api" id="toc-component-api">Component API</a>
    - <a href="#censoring" id="toc-censoring"><code>censoring</code></a>
    - <a href="#method" id="toc-method"><code>method</code></a>
    - <a href="#metric" id="toc-metric"><code>metric</code></a>

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

## Methods

Methods for assigning labels from a reference dataset to a new dataset.

| Name                                                                 | Description                    | DOI                                                  | URL                                                                                                    |
|:---------------------------------------------------------------------|:-------------------------------|:-----------------------------------------------------|:-------------------------------------------------------------------------------------------------------|
| [KNN](./methods/knn_classifier/config.vsh.yaml)                      | K-Nearest Neighbors classifier | [link](https://doi.org/10.1109/TIT.1967.1053964)     | [link](https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.KNeighborsClassifier.html)  |
| [Logistic Regression](./methods/logistic_regression/config.vsh.yaml) | Logistic regression method     |                                                      | [link](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html) |
| [Multilayer perceptron](./methods/mlp/config.vsh.yaml)               | Multilayer perceptron          | [link](https://doi.org/10.1016/0004-3702(89)90049-0) | [link](https://scikit-learn.org/stable/modules/generated/sklearn.neural_network.MLPClassifier.html)    |

## Metrics

Metrics for label projection aim to characterize how well each
classifier correctly assigns cell type labels to cells in the test set.

- **[Accuracy](./metrics/accuracy/config.vsh.yaml)**: The percentage of
  correctly predicted labels. Range: \[0, 1\]. Higher is better.
- **[F1 weighted](./metrics/f1/config.vsh.yaml)**: Calculates the F1
  score for each label, and find their average weighted by support (the
  number of true instances for each label). This alters ‘macro’ to
  account for label imbalance; it can result in an F-score that is not
  between precision and recall. Range: \[0, 1\]. Higher is better.
- **[F1 macro](./metrics/f1/config.vsh.yaml)**: Calculates the F1 score
  for each label, and find their unweighted mean. This does not take
  label imbalance into account. Range: \[0, 1\]. Higher is better.
- **[F1 micro](./metrics/f1/config.vsh.yaml)**: Calculates the F1 score
  globally by counting the total true positives, false negatives and
  false positives. Range: \[0, 1\]. Higher is better.

## Pipeline topology

``` mermaid
%%| column: screen-inset-shaded
flowchart LR
  anndata_dataset(dataset.h5ad)
  anndata_prediction(prediction.h5ad)
  anndata_score(score.h5ad)
  anndata_solution(solution.h5ad)
  anndata_test(test.h5ad)
  anndata_train(train.h5ad)
  comp_censoring[/censoring/]
  comp_method[/method/]
  comp_metric[/metric/]
  anndata_dataset---comp_censoring
  anndata_train---comp_method
  anndata_test---comp_method
  anndata_solution---comp_metric
  anndata_prediction---comp_metric
  comp_censoring-->anndata_train
  comp_censoring-->anndata_test
  comp_censoring-->anndata_solution
  comp_method-->anndata_prediction
  comp_metric-->anndata_score
```

## File format API

### `dataset.h5ad`: Preprocessed dataset

A preprocessed dataset

Used in:

- [censoring](#censoring): input (as input)

Slots:

| struct | name           | type    | description                                                         |
|:-------|:---------------|:--------|:--------------------------------------------------------------------|
| layers | counts         | integer | Raw counts                                                          |
| layers | lognorm        | double  | Log-transformed normalised counts                                   |
| obs    | label          | double  | Ground truth cell type labels                                       |
| obs    | batch          | double  | Batch information                                                   |
| uns    | dataset_id     | string  | A unique identifier for the dataset                                 |
| uns    | raw_dataset_id | string  | A unique identifier for the original dataset (before preprocessing) |

### `prediction.h5ad`: Prediction

The prediction file

Used in:

- [method](#method): output (as output)
- [metric](#metric): input_prediction (as input)

Slots:

| struct | name           | type   | description                                                         |
|:-------|:---------------|:-------|:--------------------------------------------------------------------|
| obs    | label_pred     | string | Predicted labels for the test cells.                                |
| uns    | dataset_id     | string | A unique identifier for the dataset                                 |
| uns    | raw_dataset_id | string | A unique identifier for the original dataset (before preprocessing) |
| uns    | method_id      | string | A unique identifier for the method                                  |

### `score.h5ad`: Score

Metric score file

Used in:

- [metric](#metric): output (as output)

Slots:

| struct | name           | type   | description                                                                                  |
|:-------|:---------------|:-------|:---------------------------------------------------------------------------------------------|
| uns    | dataset_id     | string | A unique identifier for the dataset                                                          |
| uns    | raw_dataset_id | string | A unique identifier for the original dataset (before preprocessing)                          |
| uns    | method_id      | string | A unique identifier for the method                                                           |
| uns    | metric_ids     | string | One or more unique metric identifiers                                                        |
| uns    | metric_values  | double | The metric values obtained for the given prediction. Must be of same length as ‘metric_ids’. |

### `solution.h5ad`: Solution

The solution for the test data

Used in:

- [censoring](#censoring): output_solution (as output)
- [metric](#metric): input_solution (as input)

Slots:

| struct | name           | type    | description                                                         |
|:-------|:---------------|:--------|:--------------------------------------------------------------------|
| layers | counts         | integer | Raw counts                                                          |
| layers | lognorm        | double  | Log-transformed normalised counts                                   |
| obs    | label          | string  | Ground truth cell type labels                                       |
| obs    | batch          | string  | Batch information                                                   |
| uns    | dataset_id     | string  | A unique identifier for the dataset                                 |
| uns    | raw_dataset_id | string  | A unique identifier for the original dataset (before preprocessing) |

### `test.h5ad`: Test data

The censored test data

Used in:

- [censoring](#censoring): output_test (as output)
- [method](#method): input_test (as input)

Slots:

| struct | name           | type    | description                                                         |
|:-------|:---------------|:--------|:--------------------------------------------------------------------|
| layers | counts         | integer | Raw counts                                                          |
| layers | lognorm        | double  | Log-transformed normalised counts                                   |
| obs    | batch          | string  | Batch information                                                   |
| uns    | dataset_id     | string  | A unique identifier for the dataset                                 |
| uns    | raw_dataset_id | string  | A unique identifier for the original dataset (before preprocessing) |

### `train.h5ad`: Training data

The training data

Used in:

- [censoring](#censoring): output_train (as output)
- [method](#method): input_train (as input)

Slots:

| struct | name           | type    | description                                                         |
|:-------|:---------------|:--------|:--------------------------------------------------------------------|
| layers | counts         | integer | Raw counts                                                          |
| layers | lognorm        | double  | Log-transformed normalised counts                                   |
| obs    | label          | string  | Ground truth cell type labels                                       |
| obs    | batch          | string  | Batch information                                                   |
| uns    | dataset_id     | string  | A unique identifier for the dataset                                 |
| uns    | raw_dataset_id | string  | A unique identifier for the original dataset (before preprocessing) |

## Component API

### `censoring`

Arguments:

| Name                | File format                     | Direction | Description          |
|:--------------------|:--------------------------------|:----------|:---------------------|
| `--input`           | [dataset.h5ad](#file-dataset)   | input     | Preprocessed dataset |
| `--output_train`    | [train.h5ad](#file-train)       | output    | Training data        |
| `--output_test`     | [test.h5ad](#file-test)         | output    | Test data            |
| `--output_solution` | [solution.h5ad](#file-solution) | output    | Solution             |

### `method`

Arguments:

| Name            | File format                         | Direction | Description   |
|:----------------|:------------------------------------|:----------|:--------------|
| `--input_train` | [train.h5ad](#file-train)           | input     | Training data |
| `--input_test`  | [test.h5ad](#file-test)             | input     | Test data     |
| `--output`      | [prediction.h5ad](#file-prediction) | output    | Prediction    |

### `metric`

Arguments:

| Name                 | File format                         | Direction | Description |
|:---------------------|:------------------------------------|:----------|:------------|
| `--input_solution`   | [solution.h5ad](#file-solution)     | input     | Solution    |
| `--input_prediction` | [prediction.h5ad](#file-prediction) | input     | Prediction  |
| `--output`           | [score.h5ad](#file-score)           | output    | Score       |
