
- <a href="#common-datasets" id="toc-common-datasets">Common datasets</a>
  - <a href="#pipeline-topology" id="toc-pipeline-topology">Pipeline
    topology</a>
  - <a href="#file-format-api" id="toc-file-format-api">File format API</a>
    - <a href="#hvg" id="toc-hvg"><code>Hvg</code></a>
    - <a href="#normalized" id="toc-normalized"><code>Normalized</code></a>
    - <a href="#pca" id="toc-pca"><code>Pca</code></a>
    - <a href="#raw" id="toc-raw"><code>Raw</code></a>
  - <a href="#component-api" id="toc-component-api">Component API</a>
    - <a href="#dataset-loader"
      id="toc-dataset-loader"><code>Dataset Loader</code></a>
    - <a href="#normalization"
      id="toc-normalization"><code>Normalization</code></a>
    - <a href="#processor-hvg"
      id="toc-processor-hvg"><code>Processor Hvg</code></a>
    - <a href="#processor-pca"
      id="toc-processor-pca"><code>Processor Pca</code></a>

# Common datasets

TODO: fill in

## Pipeline topology

``` mermaid
%%| column: screen-inset-shaded
flowchart LR
  anndata_hvg(Hvg)
  anndata_normalized(Normalized)
  anndata_pca(Pca)
  anndata_raw(Raw)
  comp_dataset_loader[/Dataset Loader/]
  comp_normalization[/Normalization/]
  comp_processor_hvg[/Processor Hvg/]
  comp_processor_pca[/Processor Pca/]
  anndata_raw---comp_normalization
  anndata_pca---comp_processor_hvg
  anndata_normalized---comp_processor_pca
  comp_dataset_loader-->anndata_raw
  comp_normalization-->anndata_normalized
  comp_processor_hvg-->anndata_hvg
  comp_processor_pca-->anndata_pca
```

## File format API

### `Hvg`

A dataset

Used in:

- [processor hvg](#processor%20hvg): output (as output)

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
| var    | hvg_ranking      | integer | A ranking of the features by hvg.                                       |
| obsm   | X_pca            | double  | The resulting PCA embedding.                                            |
| varm   | pca_loadings     | double  | The PCA loadings matrix.                                                |
| uns    | dataset_id       | string  | A unique identifier for the dataset                                     |
| uns    | normalization_id | string  | Which normalization was used                                            |
| uns    | pca_variance     | double  | The PCA variance objects.                                               |

Example:

    AnnData object
     obs: 'celltype', 'batch', 'tissue', 'size_factors'
     var: 'hvg', 'hvg_ranking'
     uns: 'dataset_id', 'normalization_id', 'pca_variance'
     obsm: 'X_pca'
     varm: 'pca_loadings'
     layers: 'counts', 'normalized'

### `Normalized`

A dataset

Used in:

- [normalization](#normalization): output (as output)
- [processor pca](#processor%20pca): input (as input)

Slots:

| struct | name             | type    | description                                                   |
|:-------|:-----------------|:--------|:--------------------------------------------------------------|
| layers | counts           | integer | Raw counts                                                    |
| layers | normalized       | double  | Normalised expression values                                  |
| obs    | celltype         | string  | Cell type information                                         |
| obs    | batch            | string  | Batch information                                             |
| obs    | tissue           | string  | Tissue information                                            |
| obs    | size_factors     | double  | The size factors created by the normalisation method, if any. |
| uns    | dataset_id       | string  | A unique identifier for the dataset                           |
| uns    | normalization_id | string  | Which normalization was used                                  |

Example:

    AnnData object
     obs: 'celltype', 'batch', 'tissue', 'size_factors'
     uns: 'dataset_id', 'normalization_id'
     layers: 'counts', 'normalized'

### `Pca`

A dataset

Used in:

- [processor hvg](#processor%20hvg): input (as input)
- [processor pca](#processor%20pca): output (as output)

Slots:

| struct | name             | type    | description                                                   |
|:-------|:-----------------|:--------|:--------------------------------------------------------------|
| layers | counts           | integer | Raw counts                                                    |
| layers | normalized       | double  | Normalised expression values                                  |
| obs    | celltype         | string  | Cell type information                                         |
| obs    | batch            | string  | Batch information                                             |
| obs    | tissue           | string  | Tissue information                                            |
| obs    | size_factors     | double  | The size factors created by the normalisation method, if any. |
| obsm   | X_pca            | double  | The resulting PCA embedding.                                  |
| varm   | pca_loadings     | double  | The PCA loadings matrix.                                      |
| uns    | dataset_id       | string  | A unique identifier for the dataset                           |
| uns    | normalization_id | string  | Which normalization was used                                  |
| uns    | pca_variance     | double  | The PCA variance objects.                                     |

Example:

    AnnData object
     obs: 'celltype', 'batch', 'tissue', 'size_factors'
     uns: 'dataset_id', 'normalization_id', 'pca_variance'
     obsm: 'X_pca'
     varm: 'pca_loadings'
     layers: 'counts', 'normalized'

### `Raw`

A raw dataset

Used in:

- [dataset loader](#dataset%20loader): output (as output)
- [normalization](#normalization): input (as input)

Slots:

| struct | name       | type    | description                         |
|:-------|:-----------|:--------|:------------------------------------|
| layers | counts     | integer | Raw counts                          |
| obs    | celltype   | string  | Cell type information               |
| obs    | batch      | string  | Batch information                   |
| obs    | tissue     | string  | Tissue information                  |
| uns    | dataset_id | string  | A unique identifier for the dataset |

Example:

    AnnData object
     obs: 'celltype', 'batch', 'tissue'
     uns: 'dataset_id'
     layers: 'counts'

## Component API

### `Dataset Loader`

Arguments:

| Name       | File format | Direction | Description |
|:-----------|:------------|:----------|:------------|
| `--output` | [Raw](#raw) | output    | Raw dataset |

### `Normalization`

Arguments:

| Name                 | File format               | Direction | Description        |
|:---------------------|:--------------------------|:----------|:-------------------|
| `--input`            | [Raw](#raw)               | input     | Raw dataset        |
| `--output`           | [Normalized](#normalized) | output    | Normalized dataset |
| `--layer_output`     | [NA](#NA)                 | input     | NA                 |
| `--obs_size_factors` | [NA](#NA)                 | input     | NA                 |

### `Processor Hvg`

Arguments:

| Name                | File format | Direction | Description     |
|:--------------------|:------------|:----------|:----------------|
| `--input`           | [Pca](#pca) | input     | Dataset+PCA     |
| `--layer_input`     | [NA](#NA)   | input     | NA              |
| `--output`          | [Hvg](#hvg) | output    | Dataset+PCA+HVG |
| `--var_hvg`         | [NA](#NA)   | input     | NA              |
| `--var_hvg_ranking` | [NA](#NA)   | input     | NA              |
| `--num_features`    | [NA](#NA)   | input     | NA              |

### `Processor Pca`

Arguments:

| Name               | File format               | Direction | Description        |
|:-------------------|:--------------------------|:----------|:-------------------|
| `--input`          | [Normalized](#normalized) | input     | Normalized dataset |
| `--layer_input`    | [NA](#NA)                 | input     | NA                 |
| `--output`         | [Pca](#pca)               | output    | Dataset+PCA        |
| `--obsm_embedding` | [NA](#NA)                 | input     | NA                 |
| `--varm_loadings`  | [NA](#NA)                 | input     | NA                 |
| `--uns_variance`   | [NA](#NA)                 | input     | NA                 |
| `--num_components` | [NA](#NA)                 | input     | NA                 |
