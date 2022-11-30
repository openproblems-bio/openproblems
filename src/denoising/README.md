
- <a href="#denoising" id="toc-denoising">Denoising</a>
  - <a href="#the-task" id="toc-the-task">The task</a>
  - <a href="#the-metrics" id="toc-the-metrics">The metrics</a>
  - <a href="#api" id="toc-api">API</a>
  - <a href="#methods" id="toc-methods">Methods</a>
  - <a href="#metrics" id="toc-metrics">Metrics</a>
  - <a href="#pipeline-topology" id="toc-pipeline-topology">Pipeline
    topology</a>
  - <a href="#file-format-api" id="toc-file-format-api">File format API</a>
    - <a href="#dataset" id="toc-dataset"><code>Dataset</code></a>
    - <a href="#denoised" id="toc-denoised"><code>Denoised</code></a>
    - <a href="#score" id="toc-score"><code>Score</code></a>
    - <a href="#test" id="toc-test"><code>Test</code></a>
    - <a href="#train" id="toc-train"><code>Train</code></a>
  - <a href="#component-api" id="toc-component-api">Component API</a>
    - <a href="#control-method"
      id="toc-control-method"><code>Control Method</code></a>
    - <a href="#method" id="toc-method"><code>Method</code></a>
    - <a href="#metric" id="toc-metric"><code>Metric</code></a>
    - <a href="#split-dataset"
      id="toc-split-dataset"><code>Split Dataset</code></a>

# Denoising

## The task

Single-cell RNA-Seq protocols only detect a fraction of the mRNA
molecules present in each cell. As a result, the measurements (UMI
counts) observed for each gene and each cell are associated with
generally high levels of technical noise ([Grün et al.,
2014](https://www.nature.com/articles/nmeth.2930)). Denoising describes
the task of estimating the true expression level of each gene in each
cell. In the single-cell literature, this task is also referred to as
*imputation*, a term which is typically used for missing data problems
in statistics. Similar to the use of the terms “dropout”, “missing
data”, and “technical zeros”, this terminology can create confusion
about the underlying measurement process ([Sarkar and Stephens,
2020](https://www.biorxiv.org/content/10.1101/2020.04.07.030007v2)).

A key challenge in evaluating denoising methods is the general lack of a
ground truth. A recent benchmark study ([Hou et al.,
2020](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02132-x))
relied on flow-sorted datasets, mixture control experiments ([Tian et
al., 2019](https://www.nature.com/articles/s41592-019-0425-8)), and
comparisons with bulk RNA-Seq data. Since each of these approaches
suffers from specific limitations, it is difficult to combine these
different approaches into a single quantitative measure of denoising
accuracy. Here, we instead rely on an approach termed molecular
cross-validation (MCV), which was specifically developed to quantify
denoising accuracy in the absence of a ground truth ([Batson et al.,
2019](https://www.biorxiv.org/content/10.1101/786269v1)). In MCV, the
observed molecules in a given scRNA-Seq dataset are first partitioned
between a *training* and a *test* dataset. Next, a denoising method is
applied to the training dataset. Finally, denoising accuracy is measured
by comparing the result to the test dataset. The authors show that both
in theory and in practice, the measured denoising accuracy is
representative of the accuracy that would be obtained on a ground truth
dataset.

## The metrics

Metrics for data denoising aim to assess denoising accuracy by comparing
the denoised *training* set to the randomly sampled *test* set.

- **MSE**: The mean squared error between the denoised counts of the
  training dataset and the true counts of the test dataset after
  reweighting by the train/test ratio.
- **Poisson**: The Poisson log likelihood of observing the true counts
  of the test dataset given the distribution given in the denoised
  dataset.

## API

Datasets should contain the raw UMI counts in `adata.X`, subsampled to
training (`adata.obsm["train"]`) and testing (`adata.obsm["test"]`)
datasets using `openproblems.tasks.denoising.datasets.utils.split_data`.

The task-specific data loader functions should split the provided raw
UMI counts into a training and a test dataset, as described by [Batson
et al., 2019](https://www.biorxiv.org/content/10.1101/786269v1). The
training dataset should be stored in `adata.obsm['train']`, and the test
dataset should be stored in `adata.obsm['test']`. Methods should store
the denoising result in `adata.obsm['denoised']`. Methods should not
edit `adata.obsm["train"]` or `adata.obsm["test"]`.

## Methods

Methods for assigning labels from a reference dataset to a new dataset.

| Name                                                                                                                                         | Type             | Description                                            | DOI                                                | URL                                               |
|:---------------------------------------------------------------------------------------------------------------------------------------------|:-----------------|:-------------------------------------------------------|:---------------------------------------------------|:--------------------------------------------------|
| [ALRA](/home/rcannood/workspace/openproblems/openproblems-v2/src/denoising/./methods/alra/config.vsh.yaml)                                   | method           | Adaptively-thresholded Low Rank Approximation (ALRA).  | [link](https://doi.org/10.1101/397588)             | [link](https://github.com/KlugerLab/ALRA)         |
| [DCA](/home/rcannood/workspace/openproblems/openproblems-v2/src/denoising/./methods/dca/config.vsh.yaml)                                     | method           | Deep Count Autoencoder                                 | [link](https://doi.org/10.1038/s41467-018-07931-2) | [link](https://github.com/theislab/dca)           |
| [knn_smooth](/home/rcannood/workspace/openproblems/openproblems-v2/src/denoising/./methods/knn_smoothing/config.vsh.yaml)                    | method           | iterative K-nearest neighbor smoothing                 | [link](https://doi.org/10.1101/217737)             | [link](https://github.com/yanailab/knn-smoothing) |
| [magic](/home/rcannood/workspace/openproblems/openproblems-v2/src/denoising/./methods/magic/config.vsh.yaml)                                 | method           | MAGIC: Markov affinity-based graph imputation of cells | [link](https://doi.org/10.1016/j.cell.2018.05.061) | [link](https://github.com/KrishnaswamyLab/MAGIC)  |
| [no denoising](/home/rcannood/workspace/openproblems/openproblems-v2/src/denoising/./control_methods/no_denoising/config.vsh.yaml)           | negative_control | negative control by copying train counts               |                                                    |                                                   |
| [perfect denoising](/home/rcannood/workspace/openproblems/openproblems-v2/src/denoising/./control_methods/perfect_denoising/config.vsh.yaml) | positive_control | Negative control by copying the train counts           |                                                    |                                                   |

## Metrics

Metrics for label projection aim to characterize how well each
classifier correctly assigns cell type labels to cells in the test set.

| Name                                                                                                             | Description                                                                                                                                                                  | Range       |
|:-----------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------------|
| [mse](/home/rcannood/workspace/openproblems/openproblems-v2/src/denoising/./metrics/mse/config.vsh.yaml)         | The mean squared error between the denoised counts of the training dataset and the true counts of the test dataset after reweighing by the train/test ratio Lower is better. | \[0, +inf\] |
| [poisson](/home/rcannood/workspace/openproblems/openproblems-v2/src/denoising/./metrics/poisson/config.vsh.yaml) | Poisson loss: measure the mean of the inconsistencies between predicted and target Lower is better.                                                                          | \[0, +inf\] |

## Pipeline topology

``` mermaid
%%| column: screen-inset-shaded
flowchart LR
  anndata_dataset(Dataset)
  anndata_denoised(Denoised)
  anndata_score(Score)
  anndata_test(Test)
  anndata_train(Train)
  comp_control_method[/Control Method/]
  comp_method[/Method/]
  comp_metric[/Metric/]
  comp_split_dataset[/Split Dataset/]
  anndata_train---comp_control_method
  anndata_test---comp_control_method
  anndata_train---comp_method
  anndata_test---comp_metric
  anndata_denoised---comp_metric
  anndata_dataset---comp_split_dataset
  comp_control_method-->anndata_denoised
  comp_method-->anndata_denoised
  comp_metric-->anndata_score
  comp_split_dataset-->anndata_train
  comp_split_dataset-->anndata_test
```

## File format API

### `Dataset`

A preprocessed dataset

Used in:

- [split dataset](#split%20dataset): input (as input)

Slots:

| struct | name       | type    | description                         |
|:-------|:-----------|:--------|:------------------------------------|
| layers | counts     | integer | Raw counts                          |
| uns    | dataset_id | string  | A unique identifier for the dataset |

Example:

    AnnData object
     uns: 'dataset_id'
     layers: 'counts'

### `Denoised`

The denoised data

Used in:

- [control method](#control%20method): output (as output)
- [method](#method): output (as output)
- [metric](#metric): input_denoised (as input)

Slots:

| struct | name       | type    | description                         |
|:-------|:-----------|:--------|:------------------------------------|
| layers | counts     | integer | Raw counts                          |
| layers | denoised   | integer | denoised data                       |
| obs    | n_counts   | string  | Raw counts                          |
| uns    | dataset_id | string  | A unique identifier for the dataset |
| uns    | method_id  | string  | A unique identifier for the method  |

Example:

    AnnData object
     obs: 'n_counts'
     uns: 'dataset_id', 'method_id'
     layers: 'counts', 'denoised'

### `Score`

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

Example:

    AnnData object
     uns: 'dataset_id', 'method_id', 'metric_ids', 'metric_values'

### `Test`

The test data

Used in:

- [control method](#control%20method): input_test (as input)
- [metric](#metric): input_test (as input)
- [split dataset](#split%20dataset): output_test (as output)

Slots:

| struct | name       | type    | description                         |
|:-------|:-----------|:--------|:------------------------------------|
| layers | counts     | integer | Raw counts                          |
| obs    | n_counts   | string  | Raw counts                          |
| uns    | dataset_id | string  | A unique identifier for the dataset |

Example:

    AnnData object
     obs: 'n_counts'
     uns: 'dataset_id'
     layers: 'counts'

### `Train`

The training data

Used in:

- [control method](#control%20method): input_train (as input)
- [method](#method): input_train (as input)
- [split dataset](#split%20dataset): output_train (as output)

Slots:

| struct | name       | type    | description                         |
|:-------|:-----------|:--------|:------------------------------------|
| layers | counts     | integer | Raw counts                          |
| obs    | n_counts   | string  | Raw counts                          |
| uns    | dataset_id | string  | A unique identifier for the dataset |

Example:

    AnnData object
     obs: 'n_counts'
     uns: 'dataset_id'
     layers: 'counts'

## Component API

### `Control Method`

Arguments:

| Name            | File format           | Direction | Description   |
|:----------------|:----------------------|:----------|:--------------|
| `--input_train` | [Train](#train)       | input     | Training data |
| `--input_test`  | [Test](#test)         | input     | Test data     |
| `--output`      | [Denoised](#denoised) | output    | Denoised data |

### `Method`

Arguments:

| Name            | File format           | Direction | Description   |
|:----------------|:----------------------|:----------|:--------------|
| `--input_train` | [Train](#train)       | input     | Training data |
| `--output`      | [Denoised](#denoised) | output    | Denoised data |

### `Metric`

Arguments:

| Name               | File format           | Direction | Description   |
|:-------------------|:----------------------|:----------|:--------------|
| `--input_test`     | [Test](#test)         | input     | Test data     |
| `--input_denoised` | [Denoised](#denoised) | input     | Denoised data |
| `--output`         | [Score](#score)       | output    | Score         |

### `Split Dataset`

Arguments:

| Name             | File format         | Direction | Description          |
|:-----------------|:--------------------|:----------|:---------------------|
| `--input`        | [Dataset](#dataset) | input     | Preprocessed dataset |
| `--output_train` | [Train](#train)     | output    | Training data        |
| `--output_test`  | [Test](#test)       | output    | Test data            |
