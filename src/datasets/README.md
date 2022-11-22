
- <a href="#common-datasets" id="toc-common-datasets">Common datasets</a>
  - <a href="#pipeline-topology" id="toc-pipeline-topology">Pipeline
    topology</a>
  - <a href="#file-format-api" id="toc-file-format-api">File format API</a>
    - <a href="#datasetpcahvg"
      id="toc-datasetpcahvg"><code>Dataset+Pca+Hvg</code></a>
    - <a href="#normalized-dataset"
      id="toc-normalized-dataset"><code>Normalized Dataset</code></a>
    - <a href="#datasetpca" id="toc-datasetpca"><code>Dataset+Pca</code></a>
    - <a href="#raw-dataset" id="toc-raw-dataset"><code>Raw Dataset</code></a>
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

## Pipeline topology

``` mermaid
%%| column: screen-inset-shaded
flowchart LR
  anndata_hvg(Dataset+Pca+Hvg)
  anndata_normalized(Normalized Dataset)
  anndata_pca(Dataset+Pca)
  anndata_raw(Raw Dataset)
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

### `Dataset+Pca+Hvg`

A normalised data with a PCA embedding and HVG selection

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
| var    | hvg_score      | integer | A ranking of the features by hvg.                                       |
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

### `Normalized Dataset`

A normalized dataset

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

### `Dataset+Pca`

A normalised data with a PCA embedding

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

### `Raw Dataset`

An unprocessed dataset as output by a dataset loader.

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

| Name       | Type                          | Direction | Description                                           |
|:-----------|:------------------------------|:----------|:------------------------------------------------------|
| `--output` | [Raw Dataset](#Raw%20dataset) | output    | An unprocessed dataset as output by a dataset loader. |

### `Normalization`

Arguments:

| Name                 | Type                                        | Direction | Description                                                  |
|:---------------------|:--------------------------------------------|:----------|:-------------------------------------------------------------|
| `--input`            | [Raw Dataset](#Raw%20dataset)               | input     | An unprocessed dataset as output by a dataset loader.        |
| `--output`           | [Normalized Dataset](#Normalized%20dataset) | output    | A normalized dataset                                         |
| `--layer_output`     | `string`                                    | input     | The name of the layer in which to store the normalized data. |
| `--obs_size_factors` | `string`                                    | input     | In which .obs slot to store the size factors (if any).       |

### `Processor Hvg`

Arguments:

| Name                | Type                                | Direction | Description                                                                |
|:--------------------|:------------------------------------|:----------|:---------------------------------------------------------------------------|
| `--input`           | [Dataset+Pca](#Dataset+PCA)         | input     | A normalised data with a PCA embedding                                     |
| `--layer_input`     | `string`                            | input     | Which layer to use as input for the PCA.                                   |
| `--output`          | [Dataset+Pca+Hvg](#Dataset+PCA+HVG) | output    | A normalised data with a PCA embedding and HVG selection                   |
| `--var_hvg`         | `string`                            | input     | In which .var slot to store whether a feature is considered to be hvg.     |
| `--var_hvg_score` | `string`                            | input     | In which .var slot to store whether a ranking of the features by variance. |
| `--num_features`    | `integer`                           | input     | The number of HVG to select                                                |

### `Processor Pca`

Arguments:

| Name               | Type                                        | Direction | Description                                                                                                          |
|:-------------------|:--------------------------------------------|:----------|:---------------------------------------------------------------------------------------------------------------------|
| `--input`          | [Normalized Dataset](#Normalized%20dataset) | input     | A normalized dataset                                                                                                 |
| `--layer_input`    | `string`                                    | input     | Which layer to use as input for the PCA.                                                                             |
| `--output`         | [Dataset+Pca](#Dataset+PCA)                 | output    | A normalised data with a PCA embedding                                                                               |
| `--obsm_embedding` | `string`                                    | input     | In which .obsm slot to store the resulting embedding.                                                                |
| `--varm_loadings`  | `string`                                    | input     | In which .varm slot to store the resulting loadings matrix.                                                          |
| `--uns_variance`   | `string`                                    | input     | In which .uns slot to store the resulting variance objects.                                                          |
| `--num_components` | `integer`                                   | input     | Number of principal components to compute. Defaults to 50, or 1 - minimum dimension size of selected representation. |
