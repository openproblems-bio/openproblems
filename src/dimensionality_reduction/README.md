
- <a href="#dimensionality-reduction"
  id="toc-dimensionality-reduction">Dimensionality reduction</a>
  - <a href="#methods" id="toc-methods">Methods</a>
  - <a href="#metrics" id="toc-metrics">Metrics</a>
  - <a href="#pipeline-topology" id="toc-pipeline-topology">Pipeline
    topology</a>
  - <a href="#file-format-api" id="toc-file-format-api">File format API</a>
    - <a href="#dataset" id="toc-dataset"><code>Dataset</code></a>
    - <a href="#reduced" id="toc-reduced"><code>Reduced</code></a>
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

# Dimensionality reduction

Dimensionality reduction is one of the key challenges in single-cell
data representation. Routine single-cell RNA sequencing (scRNA-seq)
experiments measure cells in roughly 20,000-30,000 dimensions (i.e.,
features - mostly gene transcripts but also other functional elements
encoded in mRNA such as lncRNAs). Since its inception,scRNA-seq
experiments have been growing in terms of the number of cells measured.
Originally, cutting-edge SmartSeq experiments would yield a few hundred
cells, at best. Now, it is not uncommon to see experiments that yield
over [100,000 cells](https://www.nature.com/articles/s41586-018-0590-4)
or even [\> 1 million cells.](https://doi.org/10.1126/science.aba7721).

Each *feature* in a dataset functions as a single dimension. While each
of the \~30,000 dimensions measured in each cell contribute to an
underlying data structure, the overall structure of the data is
challenging to display in few dimensions due to data sparsity and the
[*“curse of
dimensionality”*](https://en.wikipedia.org/wiki/Curse_of_dimensionality)
(distances in high dimensional data don’t distinguish data points well).
Thus, we need to find a way to [dimensionally
reduce](https://en.wikipedia.org/wiki/Dimensionality_reduction) the data
for visualization and interpretation.

## Methods

Methods to assign dimensionally-reduced 2D embedding coordinates to
adata.obsm\[‘X_emb’\].

    Warning: Unknown or uninitialised column: `paper_doi`.

    Warning: Unknown or uninitialised column: `code_url`.

| Name                                                                                                                                | Type             | Description                                                                                                                             | DOI | URL |
|:------------------------------------------------------------------------------------------------------------------------------------|:-----------------|:----------------------------------------------------------------------------------------------------------------------------------------|:----|:----|
| [densMAP](/home/jacorvar/DI/openproblems-v2/src/dimensionality_reduction/./methods/densmap/config.vsh.yaml)                         | method           | density-preserving based on UMAP                                                                                                        |     |     |
| [NeuralEE](/home/jacorvar/DI/openproblems-v2/src/dimensionality_reduction/./methods/neuralee/config.vsh.yaml)                       | method           | A neural network implementation of elastic embedding implemented in the [NeuralEE package](https://neuralee.readthedocs.io/en/latest/). |     |     |
| [PCA](/home/jacorvar/DI/openproblems-v2/src/dimensionality_reduction/./methods/pca/config.vsh.yaml)                                 | method           | Principal component analysis                                                                                                            |     |     |
| [PHATE](/home/jacorvar/DI/openproblems-v2/src/dimensionality_reduction/./methods/phate/config.vsh.yaml)                             | method           | Potential of heat-diffusion for affinity-based transition embedding                                                                     |     |     |
| [t-SNE](/home/jacorvar/DI/openproblems-v2/src/dimensionality_reduction/./methods/tsne/config.vsh.yaml)                              | method           | t-distributed stochastic neighbor embedding                                                                                             |     |     |
| [UMAP](/home/jacorvar/DI/openproblems-v2/src/dimensionality_reduction/./methods/umap/config.vsh.yaml)                               | method           | Uniform manifold approximation and projection                                                                                           |     |     |
| [Random features](/home/jacorvar/DI/openproblems-v2/src/dimensionality_reduction/./control_methods/random_features/config.vsh.yaml) | negative_control | Negative control method which generates a random embedding                                                                              |     |     |
| [True Features](/home/jacorvar/DI/openproblems-v2/src/dimensionality_reduction/./control_methods/true_features/config.vsh.yaml)     | positive_control | Positive control method which generates high-dimensional (full data) embedding                                                          |     |     |

## Metrics

Metrics for dimensionality reduction aim to compare the dimensionality
reduced dataset (the embedding) with a whole or a higher dimensional
dataset. The more similar they are, the better the reduction is.

| Name                                                                                                                        | Description                                                                                                                 | Range      |
|:----------------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------|:-----------|
| [RMSE](/home/jacorvar/DI/openproblems-v2/src/dimensionality_reduction/./metrics/rmse/config.vsh.yaml)                       | The root mean squared error between the full (or processed) data matrix and a list of dimensionally-reduced matrices NA     | \[0, NA\]  |
| [NN_Ranking](/home/jacorvar/DI/openproblems-v2/src/dimensionality_reduction/./metrics/nn_ranking/config.vsh.yaml)           | A set of metrics from the pyDRMetrics package. NA                                                                           | \[NA, NA\] |
| [Density](/home/jacorvar/DI/openproblems-v2/src/dimensionality_reduction/./metrics/density/config.vsh.yaml)                 | density preservation: correlation of local radius with the local radii in the original data space Higher is better.         | \[0, -1\]  |
| [Trustworthiness](/home/jacorvar/DI/openproblems-v2/src/dimensionality_reduction/./metrics/trustworthiness/config.vsh.yaml) | To what extent the local structure is retained in a low-dimensional embedding in a value between 0 and 1. Higher is better. | \[0, 1\]   |

## Pipeline topology

``` mermaid
%%| column: screen-inset-shaded
flowchart LR
  anndata_dataset(Dataset)
  anndata_reduced(Reduced)
  anndata_score(Score)
  anndata_test(Test)
  anndata_train(Train)
  comp_control_method[/Control Method/]
  comp_method[/Method/]
  comp_metric[/Metric/]
  comp_split_dataset[/Split Dataset/]
  anndata_dataset---comp_control_method
  anndata_train---comp_method
  anndata_reduced---comp_metric
  anndata_test---comp_metric
  anndata_dataset---comp_split_dataset
  comp_control_method-->anndata_reduced
  comp_method-->anndata_reduced
  comp_metric-->anndata_score
  comp_split_dataset-->anndata_train
  comp_split_dataset-->anndata_test
```

## File format API

### `Dataset`

A normalized data with a PCA embedding and HVG selection

Used in:

- [control method](#control%20method): input (as input)
- [split dataset](#split%20dataset): input (as input)

Slots:

| struct | name             | type    | description                                                                          |
|:-------|:-----------------|:--------|:-------------------------------------------------------------------------------------|
| layers | counts           | integer | Raw counts                                                                           |
| layers | normalized       | double  | Normalized expression values                                                         |
| obs    | celltype         | string  | Cell type information                                                                |
| obs    | batch            | string  | Batch information                                                                    |
| obs    | tissue           | string  | Tissue information                                                                   |
| obs    | size_factors     | double  | The size factors created by the normalization method, if any.                        |
| var    | hvg              | boolean | Whether or not the feature is considered to be a ‘highly variable gene’              |
| var    | hvg_score        | double  | High variability gene score (normalized dispersion). The greater, the more variable. |
| obsm   | X_pca            | double  | The resulting PCA embedding.                                                         |
| varm   | pca_loadings     | double  | The PCA loadings matrix.                                                             |
| uns    | dataset_id       | string  | A unique identifier for the dataset                                                  |
| uns    | normalization_id | string  | Which normalization was used                                                         |
| uns    | pca_variance     | double  | The PCA variance objects.                                                            |

Example:

    AnnData object
     obs: 'celltype', 'batch', 'tissue', 'size_factors'
     var: 'hvg', 'hvg_score'
     uns: 'dataset_id', 'normalization_id', 'pca_variance'
     obsm: 'X_pca'
     varm: 'pca_loadings'
     layers: 'counts', 'normalized'

### `Reduced`

A dimensionally reduced dataset

Used in:

- [control method](#control%20method): output (as output)
- [method](#method): output (as output)
- [metric](#metric): input_reduced (as input)

Slots:

| struct | name             | type   | description                          |
|:-------|:-----------------|:-------|:-------------------------------------|
| obsm   | X_emb            | double | The dimensionally reduced embedding. |
| uns    | dataset_id       | string | A unique identifier for the dataset  |
| uns    | method_id        | string | A unique identifier for the method   |
| uns    | normalization_id | string | Which normalization was used         |

Example:

    AnnData object
     uns: 'dataset_id', 'method_id', 'normalization_id'
     obsm: 'X_emb'

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

### `Test`

The test data

Used in:

- [metric](#metric): input_test (as input)
- [split dataset](#split%20dataset): output_test (as output)

Slots:

| struct | name       | type    | description                                                                          |
|:-------|:-----------|:--------|:-------------------------------------------------------------------------------------|
| layers | counts     | integer | Raw counts                                                                           |
| layers | normalized | double  | Normalized expression values                                                         |
| var    | hvg_score  | double  | High variability gene score (normalized dispersion). The greater, the more variable. |
| uns    | dataset_id | string  | A unique identifier for the dataset                                                  |

Example:

    AnnData object
     var: 'hvg_score'
     uns: 'dataset_id'
     layers: 'counts', 'normalized'

### `Train`

The training data

Used in:

- [method](#method): input (as input)
- [split dataset](#split%20dataset): output_train (as output)

Slots:

| struct | name       | type    | description                                                                          |
|:-------|:-----------|:--------|:-------------------------------------------------------------------------------------|
| layers | counts     | integer | Raw counts                                                                           |
| layers | normalized | double  | Normalized expression values                                                         |
| var    | hvg_score  | double  | High variability gene score (normalized dispersion). The greater, the more variable. |
| uns    | dataset_id | string  | A unique identifier for the dataset                                                  |

Example:

    AnnData object
     var: 'hvg_score'
     uns: 'dataset_id'
     layers: 'counts', 'normalized'

## Component API

### `Control Method`

Arguments:

| Name       | File format         | Direction | Description     |
|:-----------|:--------------------|:----------|:----------------|
| `--input`  | [Dataset](#dataset) | input     | Dataset+PCA+HVG |
| `--output` | [Reduced](#reduced) | output    | Training data   |

### `Method`

Arguments:

| Name       | File format         | Direction | Description   |
|:-----------|:--------------------|:----------|:--------------|
| `--input`  | [Train](#train)     | input     | Training data |
| `--output` | [Reduced](#reduced) | output    | Training data |

### `Metric`

Arguments:

| Name              | File format         | Direction | Description   |
|:------------------|:--------------------|:----------|:--------------|
| `--input_reduced` | [Reduced](#reduced) | input     | Training data |
| `--input_test`    | [Test](#test)       | input     | Test data     |
| `--output`        | [Score](#score)     | output    | Score         |

### `Split Dataset`

Arguments:

| Name             | File format         | Direction | Description     |
|:-----------------|:--------------------|:----------|:----------------|
| `--input`        | [Dataset](#dataset) | input     | Dataset+PCA+HVG |
| `--output_train` | [Train](#train)     | output    | Training data   |
| `--output_test`  | [Test](#test)       | output    | Test data       |
