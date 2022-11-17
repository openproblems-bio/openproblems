
- <a href="#label-projection" id="toc-label-projection">Label
  Projection</a>
  - <a href="#task-description" id="toc-task-description">Task
    description</a>
  - <a href="#methods" id="toc-methods">Methods</a>
  - <a href="#metrics" id="toc-metrics">Metrics</a>
  - <a href="#pipeline-topology" id="toc-pipeline-topology">Pipeline
    topology</a>
  - <a href="#file-format-api" id="toc-file-format-api">File format API</a>
    - <a href="#dataset.h5ad-preprocessed-dataset"
      id="toc-dataset.h5ad-preprocessed-dataset"><code>dataset.h5ad</code>:
      Preprocessed dataset</a>
    - <a href="#prediction.h5ad-prediction"
      id="toc-prediction.h5ad-prediction"><code>prediction.h5ad</code>:
      Prediction</a>
    - <a href="#score.h5ad-score"
      id="toc-score.h5ad-score"><code>score.h5ad</code>: Score</a>
    - <a href="#solution.h5ad-solution"
      id="toc-solution.h5ad-solution"><code>solution.h5ad</code>: Solution</a>
    - <a href="#test.h5ad-test-data"
      id="toc-test.h5ad-test-data"><code>test.h5ad</code>: Test data</a>
    - <a href="#train.h5ad-training-data"
      id="toc-train.h5ad-training-data"><code>train.h5ad</code>: Training
      data</a>
  - <a href="#component-api" id="toc-component-api">Component API</a>
    - <a href="#method" id="toc-method"><code>method</code></a>
    - <a href="#metric" id="toc-metric"><code>metric</code></a>
    - <a href="#split-dataset"
      id="toc-split-dataset"><code>split dataset</code></a>

# Label Projection

## Task description

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

| Name                                                                 | Description                                                                                                 | DOI                                                  | URL                                                                                                    |
|:---------------------------------------------------------------------|:------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------|:-------------------------------------------------------------------------------------------------------|
| [KNN](./methods/knn/config.vsh.yaml)                                 | K-Nearest Neighbors classifier                                                                              | [link](https://doi.org/10.1109/TIT.1967.1053964)     | [link](https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.KNeighborsClassifier.html)  |
| [Logistic Regression](./methods/logistic_regression/config.vsh.yaml) | Logistic regression method                                                                                  |                                                      | [link](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html) |
| [Multilayer perceptron](./methods/mlp/config.vsh.yaml)               | Multilayer perceptron                                                                                       | [link](https://doi.org/10.1016/0004-3702(89)90049-0) | [link](https://scikit-learn.org/stable/modules/generated/sklearn.neural_network.MLPClassifier.html)    |
| [Scanvi](./methods/scanvi/config.vsh.yaml)                           | Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models. | [link](https://doi.org/10.1101/2020.07.16.205997)    | [link](https://github.com/YosefLab/scvi-tools)                                                         |

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
  comp_method[/method/]
  comp_metric[/metric/]
  comp_split_dataset[/split dataset/]
  anndata_train---comp_method
  anndata_test---comp_method
  anndata_solution---comp_metric
  anndata_prediction---comp_metric
  anndata_dataset---comp_split_dataset
  comp_method-->anndata_prediction
  comp_metric-->anndata_score
  comp_split_dataset-->anndata_train
  comp_split_dataset-->anndata_test
  comp_split_dataset-->anndata_solution
```

## File format API

### `dataset.h5ad`: Preprocessed dataset

A preprocessed dataset

Used in:

- [split dataset](#split%20dataset): input (as input)

Slots:

| struct | name              | type    | description                                      |
|:-------|:------------------|:--------|:-------------------------------------------------|
| layers | counts            | integer | Raw counts                                       |
| layers | log_cpm           | double  | CPM normalized counts, log transformed           |
| layers | log_scran_pooling | double  | Scran pooling normalized counts, log transformed |
| layers | sqrt_cpm          | double  | CPM normalized counts, sqrt transformed          |
| obs    | label             | double  | Ground truth cell type labels                    |
| obs    | batch             | double  | Batch information                                |
| uns    | dataset_id        | string  | A unique identifier for the dataset              |

### `prediction.h5ad`: Prediction

The prediction file

Used in:

- [method](#method): output (as output)
- [metric](#metric): input_prediction (as input)

Slots:

| struct | name       | type   | description                          |
|:-------|:-----------|:-------|:-------------------------------------|
| obs    | label_pred | string | Predicted labels for the test cells. |
| uns    | dataset_id | string | A unique identifier for the dataset  |
| uns    | method_id  | string | A unique identifier for the method   |

### `score.h5ad`: Score

Metric score file

Used in:

- [metric](#metric): output (as output)

Slots:

| struct | name          | type   | description                                                                                  |
|:-------|:--------------|:-------|:---------------------------------------------------------------------------------------------|
| uns    | dataset_id    | string | A unique identifier for the dataset                                                          |
| uns    | method_id     | string | A unique identifier for the method                                                           |
| uns    | metric_ids    | string | One or more unique metric identifiers                                                        |
| uns    | metric_values | double | The metric values obtained for the given prediction. Must be of same length as ‘metric_ids’. |

### `solution.h5ad`: Solution

The solution for the test data

Used in:

- [metric](#metric): input_solution (as input)
- [split dataset](#split%20dataset): output_solution (as output)

Slots:

| struct | name              | type    | description                                      |
|:-------|:------------------|:--------|:-------------------------------------------------|
| layers | counts            | integer | Raw counts                                       |
| layers | log_cpm           | double  | CPM normalized counts, log transformed           |
| layers | log_scran_pooling | double  | Scran pooling normalized counts, log transformed |
| layers | sqrt_cpm          | double  | CPM normalized counts, sqrt transformed          |
| obs    | label             | string  | Ground truth cell type labels                    |
| obs    | batch             | string  | Batch information                                |
| uns    | dataset_id        | string  | A unique identifier for the dataset              |

### `test.h5ad`: Test data

The test data (without labels)

Used in:

- [method](#method): input_test (as input)
- [split dataset](#split%20dataset): output_test (as output)

Slots:

| struct | name              | type    | description                                      |
|:-------|:------------------|:--------|:-------------------------------------------------|
| layers | counts            | integer | Raw counts                                       |
| layers | log_cpm           | double  | CPM normalized counts, log transformed           |
| layers | log_scran_pooling | double  | Scran pooling normalized counts, log transformed |
| layers | sqrt_cpm          | double  | CPM normalized counts, sqrt transformed          |
| obs    | batch             | string  | Batch information                                |
| uns    | dataset_id        | string  | A unique identifier for the dataset              |

### `train.h5ad`: Training data

The training data

Used in:

- [method](#method): input_train (as input)
- [split dataset](#split%20dataset): output_train (as output)

Slots:

| struct | name              | type    | description                                      |
|:-------|:------------------|:--------|:-------------------------------------------------|
| layers | counts            | integer | Raw counts                                       |
| layers | log_cpm           | double  | CPM normalized counts, log transformed           |
| layers | log_scran_pooling | double  | Scran pooling normalized counts, log transformed |
| layers | sqrt_cpm          | double  | CPM normalized counts, sqrt transformed          |
| obs    | label             | string  | Ground truth cell type labels                    |
| obs    | batch             | string  | Batch information                                |
| uns    | dataset_id        | string  | A unique identifier for the dataset              |

## Component API

### `method`

Arguments:

| Name            | File format     | Direction | Description   |
|:----------------|:----------------|:----------|:--------------|
| `--input_train` | train.h5ad      | input     | Training data |
| `--input_test`  | test.h5ad       | input     | Test data     |
| `--output`      | prediction.h5ad | output    | Prediction    |

### `metric`

Arguments:

| Name                 | File format     | Direction | Description |
|:---------------------|:----------------|:----------|:------------|
| `--input_solution`   | solution.h5ad   | input     | Solution    |
| `--input_prediction` | prediction.h5ad | input     | Prediction  |
| `--output`           | score.h5ad      | output    | Score       |

### `split dataset`

Arguments:

| Name                | File format   | Direction | Description          |
|:--------------------|:--------------|:----------|:---------------------|
| `--input`           | dataset.h5ad  | input     | Preprocessed dataset |
| `--output_train`    | train.h5ad    | output    | Training data        |
| `--output_test`     | test.h5ad     | output    | Test data            |
| `--output_solution` | solution.h5ad | output    | Solution             |
