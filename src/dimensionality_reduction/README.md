
- <a href="#label-projection" id="toc-label-projection">Label
  Projection</a>
  - <a href="#methods" id="toc-methods">Methods</a>
  - <a href="#metrics" id="toc-metrics">Metrics</a>
  - <a href="#pipeline-topology" id="toc-pipeline-topology">Pipeline
    topology</a>
  - <a href="#file-format-api" id="toc-file-format-api">File format API</a>
    - <a href="#dataset" id="toc-dataset"><code>Dataset</code></a>
    - <a href="#reduced" id="toc-reduced"><code>Reduced</code></a>
    - <a href="#score" id="toc-score"><code>Score</code></a>
  - <a href="#component-api" id="toc-component-api">Component API</a>
    - <a href="#control-method"
      id="toc-control-method"><code>Control Method</code></a>
    - <a href="#method" id="toc-method"><code>Method</code></a>
    - <a href="#metric" id="toc-metric"><code>Metric</code></a>

# Label Projection

## Methods

Methods for assigning labels from a reference dataset to a new dataset.

    Warning: Unknown or uninitialised column: `paper_doi`.

    Warning: Unknown or uninitialised column: `code_url`.

| Name                                                                                                                                                      | Type             | Description                                                            | DOI | URL |
|:----------------------------------------------------------------------------------------------------------------------------------------------------------|:-----------------|:-----------------------------------------------------------------------|:----|:----|
| [densMAP](/home/rcannood/workspace/openproblems/openproblems-v2/src/dimensionality_reduction/./methods/densmap/config.vsh.yaml)                           | method           | density-preserving based on UMAP                                       |     |     |
| [PHATE](/home/rcannood/workspace/openproblems/openproblems-v2/src/dimensionality_reduction/./methods/phate/config.vsh.yaml)                               | method           |                                                                        |     |     |
| [t-SNE](/home/rcannood/workspace/openproblems/openproblems-v2/src/dimensionality_reduction/./methods/tsne/config.vsh.yaml)                                | method           | t-distributed stochastic neighbor embedding                            |     |     |
| [UMAP](/home/rcannood/workspace/openproblems/openproblems-v2/src/dimensionality_reduction/./methods/umap/config.vsh.yaml)                                 | method           | Uniform manifold approximation and projection                          |     |     |
| [Random features](/home/rcannood/workspace/openproblems/openproblems-v2/src/dimensionality_reduction/./control_methods/random_features/config.vsh.yaml)   | negative_control | Negative control method which generates a random embedding             |     |     |
| [High-dimensional PCA](/home/rcannood/workspace/openproblems/openproblems-v2/src/dimensionality_reduction/./control_methods/high_dim_pca/config.vsh.yaml) | positive_control | Positive control method which generates high-dimensional PCA embedding |     |     |

## Metrics

Metrics for label projection aim to characterize how well each
classifier correctly assigns cell type labels to cells in the test set.

    Warning: Unknown or uninitialised column: `description`.

    Warning: Unknown or uninitialised column: `maximize`.

    Warning: Unknown or uninitialised column: `min`.

    Warning: Unknown or uninitialised column: `max`.

| Name                                                                                                                      | Description | Range      |
|:--------------------------------------------------------------------------------------------------------------------------|:------------|:-----------|
| [RMSE](/home/rcannood/workspace/openproblems/openproblems-v2/src/dimensionality_reduction/./metrics/rmse/config.vsh.yaml) | NA NA       | \[NA, NA\] |

## Pipeline topology

``` mermaid
%%| column: screen-inset-shaded
flowchart LR
  anndata_dataset(Dataset)
  anndata_reduced(Reduced)
  anndata_score(Score)
  comp_control_method[/Control Method/]
  comp_method[/Method/]
  comp_metric[/Metric/]
  anndata_dataset---comp_control_method
  anndata_dataset---comp_method
  anndata_reduced---comp_metric
  comp_control_method-->anndata_reduced
  comp_method-->anndata_reduced
  comp_metric-->anndata_score
```

## File format API

### `Dataset`

A normalised data with a PCA embedding and HVG selection

Used in:

- [control method](#control%20method): input (as input)
- [method](#method): input (as input)

Slots:

| struct | name             | type    | description                                                             |
|:-------|:-----------------|:--------|:------------------------------------------------------------------------|
| layers | counts           | integer | Raw counts                                                              |
| layers | normalized       | double  | Normalized expression values                                            |
| obs    | celltype         | string  | Cell type information                                                   |
| obs    | batch            | string  | Batch information                                                       |
| obs    | tissue           | string  | Tissue information                                                      |
| obs    | size_factors     | double  | The size factors created by the normalization method, if any.           |
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

### `Reduced`

A dimensionality reduced dataset

Used in:

- [control method](#control%20method): output (as output)
- [method](#method): output (as output)
- [metric](#metric): input (as input)

Slots:

| struct | name             | type    | description                                                             |
|:-------|:-----------------|:--------|:------------------------------------------------------------------------|
| layers | counts           | integer | Raw counts                                                              |
| layers | normalized       | double  | Normalized expression values                                            |
| obs    | celltype         | string  | Cell type information                                                   |
| obs    | batch            | string  | Batch information                                                       |
| obs    | tissue           | string  | Tissue information                                                      |
| obs    | size_factors     | double  | The size factors created by the normalization method, if any.           |
| var    | hvg              | boolean | Whether or not the feature is considered to be a ‘highly variable gene’ |
| var    | hvg_score        | integer | A ranking of the features by hvg.                                       |
| obsm   | X_pca            | double  | The resulting PCA embedding.                                            |
| obsm   | X_emb            | double  | The resulting t-SNE embedding.                                          |
| uns    | dataset_id       | string  | A unique identifier for the dataset                                     |
| uns    | method_id        | string  | A unique identifier for the method                                      |
| uns    | normalization_id | string  | Which normalization was used                                            |

Example:

    AnnData object
     obs: 'celltype', 'batch', 'tissue', 'size_factors'
     var: 'hvg', 'hvg_score'
     uns: 'dataset_id', 'method_id', 'normalization_id'
     obsm: 'X_pca', 'X_emb'
     layers: 'counts', 'normalized'

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

## Component API

### `Control Method`

Arguments:

| Name       | File format         | Direction | Description |
|:-----------|:--------------------|:----------|:------------|
| `--input`  | [Dataset](#dataset) | input     | NA          |
| `--output` | [Reduced](#reduced) | output    | NA          |

### `Method`

Arguments:

| Name       | File format         | Direction | Description |
|:-----------|:--------------------|:----------|:------------|
| `--input`  | [Dataset](#dataset) | input     | NA          |
| `--output` | [Reduced](#reduced) | output    | NA          |

### `Metric`

Arguments:

| Name       | File format         | Direction | Description |
|:-----------|:--------------------|:----------|:------------|
| `--input`  | [Reduced](#reduced) | input     | NA          |
| `--output` | [Score](#score)     | output    | Score       |
