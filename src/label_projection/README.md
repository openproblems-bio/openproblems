
- <a href="#label-projection" id="toc-label-projection">Label
  Projection</a>
  - <a href="#task-description" id="toc-task-description">Task
    description</a>
  - <a href="#methods" id="toc-methods">Methods</a>
  - <a href="#metrics" id="toc-metrics">Metrics</a>
  - <a href="#pipeline-topology" id="toc-pipeline-topology">Pipeline
    topology</a>
  - <a href="#file-format-api" id="toc-file-format-api">File format API</a>
    - <a href="#dataset" id="toc-dataset"><code>Dataset</code></a>
    - <a href="#prediction" id="toc-prediction"><code>Prediction</code></a>
    - <a href="#score" id="toc-score"><code>Score</code></a>
    - <a href="#solution" id="toc-solution"><code>Solution</code></a>
    - <a href="#test" id="toc-test"><code>Test</code></a>
    - <a href="#train" id="toc-train"><code>Train</code></a>
  - <a href="#component-api" id="toc-component-api">Component API</a>
    - <a href="#control-method"
      id="toc-control-method"><code>Control Method</code></a>
    - <a href="#method" id="toc-method"><code>Method</code></a>
    - <a href="#metric" id="toc-metric"><code>Metric</code></a>
    - <a href="#split-dataset"
      id="toc-split-dataset"><code>Split Dataset</code></a>

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

| Name                                                                 | Type             | Description                                           | DOI                                                  | URL                                                                                                    |
|:---------------------------------------------------------------------|:-----------------|:------------------------------------------------------|:-----------------------------------------------------|:-------------------------------------------------------------------------------------------------------|
| [KNN](./methods/knn/config.vsh.yaml)                                 | method           | K-Nearest Neighbors classifier                        | [link](https://doi.org/10.1109/TIT.1967.1053964)     | [link](https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.KNeighborsClassifier.html)  |
| [Logistic Regression](./methods/logistic_regression/config.vsh.yaml) | method           | Logistic regression method                            |                                                      | [link](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html) |
| [Multilayer perceptron](./methods/mlp/config.vsh.yaml)               | method           | Multilayer perceptron                                 | [link](https://doi.org/10.1016/0004-3702(89)90049-0) | [link](https://scikit-learn.org/stable/modules/generated/sklearn.neural_network.MLPClassifier.html)    |
| [Majority Vote](./control_methods/majority_vote/config.vsh.yaml)     | negative_control | Baseline method using majority voting                 |                                                      |                                                                                                        |
| [Random Labels](./control_methods/random_labels/config.vsh.yaml)     | negative_control | Negative control method which generates random labels |                                                      |                                                                                                        |
| [True labels](./control_methods/true_labels/config.vsh.yaml)         | positive_control | Positive control method by returning the true labels  |                                                      |                                                                                                        |

## Metrics

Metrics for label projection aim to characterize how well each
classifier correctly assigns cell type labels to cells in the test set.

| Name                                           | Description                                                                                                                                                                                                                                                                   | Range    |
|:-----------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:---------|
| [Accuracy](./metrics/accuracy/config.vsh.yaml) | The percentage of correctly predicted labels. Higher is better.                                                                                                                                                                                                               | \[0, 1\] |
| [F1 weighted](./metrics/f1/config.vsh.yaml)    | Calculates the F1 score for each label, and find their average weighted by support (the number of true instances for each label). This alters ‘macro’ to account for label imbalance; it can result in an F-score that is not between precision and recall. Higher is better. | \[0, 1\] |
| [F1 macro](./metrics/f1/config.vsh.yaml)       | Calculates the F1 score for each label, and find their unweighted mean. This does not take label imbalance into account. Higher is better.                                                                                                                                    | \[0, 1\] |
| [F1 micro](./metrics/f1/config.vsh.yaml)       | Calculates the F1 score globally by counting the total true positives, false negatives and false positives. Higher is better.                                                                                                                                                 | \[0, 1\] |

## Pipeline topology

``` mermaid
%%| column: screen-inset-shaded
flowchart LR
  anndata_dataset(Dataset)
  anndata_prediction(Prediction)
  anndata_score(Score)
  anndata_solution(Solution)
  anndata_test(Test)
  anndata_train(Train)
  comp_control_method[/Control Method/]
  comp_method[/Method/]
  comp_metric[/Metric/]
  comp_split_dataset[/Split Dataset/]
  anndata_train---comp_control_method
  anndata_test---comp_control_method
  anndata_solution---comp_control_method
  anndata_train---comp_method
  anndata_test---comp_method
  anndata_solution---comp_metric
  anndata_prediction---comp_metric
  anndata_dataset---comp_split_dataset
  comp_control_method-->anndata_prediction
  comp_method-->anndata_prediction
  comp_metric-->anndata_score
  comp_split_dataset-->anndata_train
  comp_split_dataset-->anndata_test
  comp_split_dataset-->anndata_solution
```

## File format API

### `Dataset`

A normalised data with a PCA embedding and HVG selection

Used in:

- [split dataset](#split%20dataset): input (as input)

Slots:

| struct | name             | type    | description                                                             |
|:-------|:-----------------|:--------|:------------------------------------------------------------------------|
| layers | counts           | integer | Raw counts                                                              |
| layers | normalized       | double  | Normalised expression values                                            |
| obs    | celltype         | string  | Cell type information                                                   |
| obs    | batch            | string  | Batch information                                                       |
| obs    | tissue           | string  | Tissue information                                                      |
| obs    | size_factors     | double  | The size factors created by the normalisation method, if any.           |
| var    | hvg              | boolean | Whether or not the feature is considered to be a ‘highly variable gene’ |
| var    | hvg_score        | integer | A ranking of the features by hvg.                                       |
| obsm   | X_pca            | double  | The resulting PCA embedding.                                            |
| varm   | pca_loadings     | double  | The PCA loadings matrix.                                                |
| uns    | dataset_id       | string  | A unique identifier for the dataset                                     |
| uns    | normalization_id | string  | Which normalization was used                                            |
| uns    | pca_variance     | double  | The PCA variance objects.                                               |

Example:

    AnnData object
     obs: 'celltype', 'batch', 'tissue', 'size_factors'
     var: 'hvg', 'hvg_score'
     uns: 'dataset_id', 'normalization_id', 'pca_variance'
     obsm: 'X_pca'
     varm: 'pca_loadings'
     layers: 'counts', 'normalized'

### `Prediction`

The prediction file

Used in:

- [control method](#control%20method): output (as output)
- [method](#method): output (as output)
- [metric](#metric): input_prediction (as input)

Slots:

| struct | name             | type   | description                          |
|:-------|:-----------------|:-------|:-------------------------------------|
| obs    | label_pred       | string | Predicted labels for the test cells. |
| uns    | dataset_id       | string | A unique identifier for the dataset  |
| uns    | normalization_id | string | Which normalization was used         |
| uns    | method_id        | string | A unique identifier for the method   |

Example:

    AnnData object
     obs: 'label_pred'
     uns: 'dataset_id', 'normalization_id', 'method_id'

### `Score`

Metric score file

Used in:

- [metric](#metric): output (as output)

Slots:

| struct | name             | type   | description                                                                                  |
|:-------|:-----------------|:-------|:---------------------------------------------------------------------------------------------|
| uns    | dataset_id       | string | A unique identifier for the dataset                                                          |
| uns    | normalization_id | string | Which normalization was used                                                                 |
| uns    | method_id        | string | A unique identifier for the method                                                           |
| uns    | metric_ids       | string | One or more unique metric identifiers                                                        |
| uns    | metric_values    | double | The metric values obtained for the given prediction. Must be of same length as ‘metric_ids’. |

Example:

    AnnData object
     uns: 'dataset_id', 'normalization_id', 'method_id', 'metric_ids', 'metric_values'

### `Solution`

The solution for the test data

Used in:

- [control method](#control%20method): input_solution (as input)
- [metric](#metric): input_solution (as input)
- [split dataset](#split%20dataset): output_solution (as output)

Slots:

| struct | name             | type    | description                                                             |
|:-------|:-----------------|:--------|:------------------------------------------------------------------------|
| layers | counts           | integer | Raw counts                                                              |
| layers | normalized       | double  | Normalized counts                                                       |
| obs    | label            | string  | Ground truth cell type labels                                           |
| obs    | batch            | string  | Batch information                                                       |
| var    | hvg              | boolean | Whether or not the feature is considered to be a ‘highly variable gene’ |
| var    | hvg_score        | integer | A ranking of the features by hvg.                                       |
| obsm   | X_pca            | double  | The resulting PCA embedding.                                            |
| uns    | dataset_id       | string  | A unique identifier for the dataset                                     |
| uns    | normalization_id | string  | Which normalization was used                                            |

Example:

    AnnData object
     obs: 'label', 'batch'
     var: 'hvg', 'hvg_score'
     uns: 'dataset_id', 'normalization_id'
     obsm: 'X_pca'
     layers: 'counts', 'normalized'

### `Test`

The test data (without labels)

Used in:

- [control method](#control%20method): input_test (as input)
- [method](#method): input_test (as input)
- [split dataset](#split%20dataset): output_test (as output)

Slots:

| struct | name             | type    | description                                                             |
|:-------|:-----------------|:--------|:------------------------------------------------------------------------|
| layers | counts           | integer | Raw counts                                                              |
| layers | normalized       | double  | Normalized counts                                                       |
| obs    | batch            | string  | Batch information                                                       |
| var    | hvg              | boolean | Whether or not the feature is considered to be a ‘highly variable gene’ |
| var    | hvg_score        | integer | A ranking of the features by hvg.                                       |
| obsm   | X_pca            | double  | The resulting PCA embedding.                                            |
| uns    | dataset_id       | string  | A unique identifier for the dataset                                     |
| uns    | normalization_id | string  | Which normalization was used                                            |

Example:

    AnnData object
     obs: 'batch'
     var: 'hvg', 'hvg_score'
     uns: 'dataset_id', 'normalization_id'
     obsm: 'X_pca'
     layers: 'counts', 'normalized'

### `Train`

The training data

Used in:

- [control method](#control%20method): input_train (as input)
- [method](#method): input_train (as input)
- [split dataset](#split%20dataset): output_train (as output)

Slots:

| struct | name             | type    | description                                                             |
|:-------|:-----------------|:--------|:------------------------------------------------------------------------|
| layers | counts           | integer | Raw counts                                                              |
| layers | normalized       | double  | Normalized counts                                                       |
| obs    | label            | string  | Ground truth cell type labels                                           |
| obs    | batch            | string  | Batch information                                                       |
| var    | hvg              | boolean | Whether or not the feature is considered to be a ‘highly variable gene’ |
| var    | hvg_score        | integer | A ranking of the features by hvg.                                       |
| obsm   | X_pca            | double  | The resulting PCA embedding.                                            |
| uns    | dataset_id       | string  | A unique identifier for the dataset                                     |
| uns    | normalization_id | string  | Which normalization was used                                            |

Example:

    AnnData object
     obs: 'label', 'batch'
     var: 'hvg', 'hvg_score'
     uns: 'dataset_id', 'normalization_id'
     obsm: 'X_pca'
     layers: 'counts', 'normalized'

## Component API

### `Control Method`

Arguments:

| Name               | File format               | Direction | Description   |
|:-------------------|:--------------------------|:----------|:--------------|
| `--input_train`    | [Train](#train)           | input     | Training data |
| `--input_test`     | [Test](#test)             | input     | Test data     |
| `--input_solution` | [Solution](#solution)     | input     | Solution      |
| `--output`         | [Prediction](#prediction) | output    | Prediction    |

### `Method`

Arguments:

| Name            | File format               | Direction | Description   |
|:----------------|:--------------------------|:----------|:--------------|
| `--input_train` | [Train](#train)           | input     | Training data |
| `--input_test`  | [Test](#test)             | input     | Test data     |
| `--output`      | [Prediction](#prediction) | output    | Prediction    |

### `Metric`

Arguments:

| Name                 | File format               | Direction | Description |
|:---------------------|:--------------------------|:----------|:------------|
| `--input_solution`   | [Solution](#solution)     | input     | Solution    |
| `--input_prediction` | [Prediction](#prediction) | input     | Prediction  |
| `--output`           | [Score](#score)           | output    | Score       |

### `Split Dataset`

Arguments:

| Name                | File format           | Direction | Description   |
|:--------------------|:----------------------|:----------|:--------------|
| `--input`           | [Dataset](#dataset)   | input     | NA            |
| `--output_train`    | [Train](#train)       | output    | Training data |
| `--output_test`     | [Test](#test)         | output    | Test data     |
| `--output_solution` | [Solution](#solution) | output    | Solution      |
