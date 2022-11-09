
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

## API

``` mermaid
%%| column: screen-inset-shaded
flowchart LR
  anndata_dataset(dataset)
  anndata_prediction(prediction)
  anndata_score(score)
  anndata_solution(solution)
  anndata_test(test)
  anndata_train(train)
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

### dataset

Used in:

- censoring: input (as input)

Slots:

| struct | name           | type    | description                                                         |
|:-------|:---------------|:--------|:--------------------------------------------------------------------|
| layers | counts         | integer | Raw counts                                                          |
| layers | lognorm        | double  | Log-transformed normalised counts                                   |
| obs    | label          | double  | Ground truth cell type labels                                       |
| obs    | batch          | double  | Batch information                                                   |
| uns    | dataset_id     | string  | A unique identifier for the dataset                                 |
| uns    | raw_dataset_id | string  | A unique identifier for the original dataset (before preprocessing) |

### prediction

Used in:

- method: output (as output)
- metric: input_prediction (as input)

Slots:

| struct | name           | type   | description                                                         |
|:-------|:---------------|:-------|:--------------------------------------------------------------------|
| obs    | label_pred     | string | Predicted labels for the test cells.                                |
| uns    | dataset_id     | string | A unique identifier for the dataset                                 |
| uns    | raw_dataset_id | string | A unique identifier for the original dataset (before preprocessing) |
| uns    | method_id      | string | A unique identifier for the method                                  |

### score

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

### solution

Used in:

- censoring: output_solution (as output)
- metric: input_solution (as input)

Slots:

| struct | name           | type    | description                                                         |
|:-------|:---------------|:--------|:--------------------------------------------------------------------|
| layers | counts         | integer | Raw counts                                                          |
| layers | lognorm        | double  | Log-transformed normalised counts                                   |
| obs    | label          | string  | Ground truth cell type labels                                       |
| obs    | batch          | string  | Batch information                                                   |
| uns    | dataset_id     | string  | A unique identifier for the dataset                                 |
| uns    | raw_dataset_id | string  | A unique identifier for the original dataset (before preprocessing) |

### test

Used in:

- censoring: output_test (as output)
- method: input_test (as input)

Slots:

| struct | name           | type    | description                                                         |
|:-------|:---------------|:--------|:--------------------------------------------------------------------|
| layers | counts         | integer | Raw counts                                                          |
| layers | lognorm        | double  | Log-transformed normalised counts                                   |
| obs    | batch          | string  | Batch information                                                   |
| uns    | dataset_id     | string  | A unique identifier for the dataset                                 |
| uns    | raw_dataset_id | string  | A unique identifier for the original dataset (before preprocessing) |

### train

Used in:

- censoring: output_train (as output)
- method: input_train (as input)

Slots:

| struct | name           | type    | description                                                         |
|:-------|:---------------|:--------|:--------------------------------------------------------------------|
| layers | counts         | integer | Raw counts                                                          |
| layers | lognorm        | double  | Log-transformed normalised counts                                   |
| obs    | label          | string  | Ground truth cell type labels                                       |
| obs    | batch          | string  | Batch information                                                   |
| uns    | dataset_id     | string  | A unique identifier for the dataset                                 |
| uns    | raw_dataset_id | string  | A unique identifier for the original dataset (before preprocessing) |
