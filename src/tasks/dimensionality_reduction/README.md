# Dimensionality reduction for 2D visualization


Reduction of high-dimensional datasets to 2D for visualization &
interpretation

Path:
[`src/tasks/dimensionality_reduction`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/tasks/dimensionality_reduction)

## Motivation

Data visualisation is an important part of all stages of single-cell
analysis, from initial quality control to interpretation and
presentation of final results. For bulk RNA-seq studies, linear
dimensionality reduction techniques such as PCA and MDS are commonly
used to visualise the variation between samples. While these methods are
highly effective they can only be used to show the first few components
of variation which cannot fully represent the increased complexity and
number of observations in single-cell datasets. For this reason
non-linear techniques (most notably t-SNE and UMAP) have become the
standard for visualising single-cell studies. These methods attempt to
compress a dataset into a two-dimensional space while attempting to
capture as much of the variance between observations as possible. Many
methods for solving this problem now exist. In general these methods try
to preserve distances, while some additionally consider aspects such as
density within the embedded space or conservation of continuous
trajectories. Despite almost every single-cell study using one of these
visualisations there has been debate as to whether they can effectively
capture the variation in single-cell datasets \[@chari2023speciousart\].

## Description

The dimensionality reduction task attempts to quantify the ability of
methods to embed the information present in complex single-cell studies
into a two-dimensional space. Thus, this task is specifically designed
for dimensionality reduction for visualisation and does not consider
other uses of dimensionality reduction in standard single-cell workflows
such as improving the signal-to-noise ratio (and in fact several of the
methods use PCA as a pre-processing step for this reason). Unlike most
tasks, methods for the dimensionality reduction task must accept a
matrix containing expression values normalised to 10,000 counts per cell
and log transformed (log-10k) and produce a two-dimensional coordinate
for each cell. Pre-normalised matrices are required to enforce
consistency between the metric evaluation (which generally requires
normalised data) and the method runs. When these are not consistent,
methods that use the same normalisation as used in the metric tend to
score more highly. For some methods we also evaluate the pre-processing
recommended by the method.

## Authors & contributors

| name                   | roles              |
|:-----------------------|:-------------------|
| Luke Zappia            | maintainer, author |
| Michal Klein           | author             |
| Scott Gigante          | author             |
| Ben DeMeo              | author             |
| Robrecht Cannoodt      | author             |
| Kai Waldrant           | contributor        |
| Sai Nirmayi Yasa       | contributor        |
| Juan A. Cordero Varela | contributor        |

## API

``` mermaid
flowchart LR
  file_common_dataset("Common dataset")
  comp_process_dataset[/"Data processor"/]
  file_dataset("Dataset")
  file_solution("Test data")
  comp_control_method[/"Control method"/]
  comp_method[/"Method"/]
  comp_metric[/"Metric"/]
  file_embedding("Embedding")
  file_score("Score")
  file_common_dataset---comp_process_dataset
  comp_process_dataset-->file_dataset
  comp_process_dataset-->file_solution
  file_dataset---comp_control_method
  file_dataset---comp_method
  file_solution---comp_control_method
  file_solution---comp_metric
  comp_control_method-->file_embedding
  comp_method-->file_embedding
  comp_metric-->file_score
  file_embedding---comp_metric
```

## File format: Common dataset

A dataset processed by the common dataset processing pipeline.

Example file: `resources_test/common/pancreas/dataset.h5ad`

Description:

This dataset contains both raw counts and normalized data matrices, as
well as a PCA embedding, HVG selection and a kNN graph.

Format:

<div class="small">

    AnnData object
     obs: 'dataset_id', 'assay', 'assay_ontology_term_id', 'cell_type', 'cell_type_ontology_term_id', 'development_stage', 'development_stage_ontology_term_id', 'disease', 'disease_ontology_term_id', 'donor_id', 'is_primary_data', 'organism', 'organism_ontology_term_id', 'self_reported_ethnicity', 'self_reported_ethnicity_ontology_term_id', 'sex', 'sex_ontology_term_id', 'suspension_type', 'tissue', 'tissue_ontology_term_id', 'tissue_general', 'tissue_general_ontology_term_id', 'batch', 'soma_joinid', 'size_factors'
     var: 'feature_id', 'feature_name', 'soma_joinid', 'hvg', 'hvg_score'
     obsm: 'X_pca'
     obsp: 'knn_distances', 'knn_connectivities'
     varm: 'pca_loadings'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'dataset_name', 'dataset_url', 'dataset_reference', 'dataset_summary', 'dataset_description', 'dataset_organism', 'normalization_id', 'pca_variance', 'knn'

</div>

Slot description:

<div class="small">

| Slot                                              | Type      | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
|:--------------------------------------------------|:----------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `obs["dataset_id"]`                               | `string`  | (*Optional*) Identifier for the dataset from which the cell data is derived, useful for tracking and referencing purposes.                                                                                                                                                                                                                                                                                                                                                                    |
| `obs["assay"]`                                    | `string`  | (*Optional*) Type of assay used to generate the cell data, indicating the methodology or technique employed.                                                                                                                                                                                                                                                                                                                                                                                  |
| `obs["assay_ontology_term_id"]`                   | `string`  | (*Optional*) Experimental Factor Ontology (`EFO:`) term identifier for the assay, providing a standardized reference to the assay type.                                                                                                                                                                                                                                                                                                                                                       |
| `obs["cell_type"]`                                | `string`  | (*Optional*) Classification of the cell type based on its characteristics and function within the tissue or organism.                                                                                                                                                                                                                                                                                                                                                                         |
| `obs["cell_type_ontology_term_id"]`               | `string`  | (*Optional*) Cell Ontology (`CL:`) term identifier for the cell type, offering a standardized reference to the specific cell classification.                                                                                                                                                                                                                                                                                                                                                  |
| `obs["development_stage"]`                        | `string`  | (*Optional*) Stage of development of the organism or tissue from which the cell is derived, indicating its maturity or developmental phase.                                                                                                                                                                                                                                                                                                                                                   |
| `obs["development_stage_ontology_term_id"]`       | `string`  | (*Optional*) Ontology term identifier for the developmental stage, providing a standardized reference to the organism’s developmental phase. If the organism is human (`organism_ontology_term_id == 'NCBITaxon:9606'`), then the Human Developmental Stages (`HsapDv:`) ontology is used. If the organism is mouse (`organism_ontology_term_id == 'NCBITaxon:10090'`), then the Mouse Developmental Stages (`MmusDv:`) ontology is used. Otherwise, the Uberon (`UBERON:`) ontology is used. |
| `obs["disease"]`                                  | `string`  | (*Optional*) Information on any disease or pathological condition associated with the cell or donor.                                                                                                                                                                                                                                                                                                                                                                                          |
| `obs["disease_ontology_term_id"]`                 | `string`  | (*Optional*) Ontology term identifier for the disease, enabling standardized disease classification and referencing. Must be a term from the Mondo Disease Ontology (`MONDO:`) ontology term, or `PATO:0000461` from the Phenotype And Trait Ontology (`PATO:`).                                                                                                                                                                                                                              |
| `obs["donor_id"]`                                 | `string`  | (*Optional*) Identifier for the donor from whom the cell sample is obtained.                                                                                                                                                                                                                                                                                                                                                                                                                  |
| `obs["is_primary_data"]`                          | `boolean` | (*Optional*) Indicates whether the data is primary (directly obtained from experiments) or has been computationally derived from other primary data.                                                                                                                                                                                                                                                                                                                                          |
| `obs["organism"]`                                 | `string`  | (*Optional*) Organism from which the cell sample is obtained.                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `obs["organism_ontology_term_id"]`                | `string`  | (*Optional*) Ontology term identifier for the organism, providing a standardized reference for the organism. Must be a term from the NCBI Taxonomy Ontology (`NCBITaxon:`) which is a child of `NCBITaxon:33208`.                                                                                                                                                                                                                                                                             |
| `obs["self_reported_ethnicity"]`                  | `string`  | (*Optional*) Ethnicity of the donor as self-reported, relevant for studies considering genetic diversity and population-specific traits.                                                                                                                                                                                                                                                                                                                                                      |
| `obs["self_reported_ethnicity_ontology_term_id"]` | `string`  | (*Optional*) Ontology term identifier for the self-reported ethnicity, providing a standardized reference for ethnic classifications. If the organism is human (`organism_ontology_term_id == 'NCBITaxon:9606'`), then the Human Ancestry Ontology (`HANCESTRO:`) is used.                                                                                                                                                                                                                    |
| `obs["sex"]`                                      | `string`  | (*Optional*) Biological sex of the donor or source organism, crucial for studies involving sex-specific traits or conditions.                                                                                                                                                                                                                                                                                                                                                                 |
| `obs["sex_ontology_term_id"]`                     | `string`  | (*Optional*) Ontology term identifier for the biological sex, ensuring standardized classification of sex. Only `PATO:0000383`, `PATO:0000384` and `PATO:0001340` are allowed.                                                                                                                                                                                                                                                                                                                |
| `obs["suspension_type"]`                          | `string`  | (*Optional*) Type of suspension or medium in which the cells were stored or processed, important for understanding cell handling and conditions.                                                                                                                                                                                                                                                                                                                                              |
| `obs["tissue"]`                                   | `string`  | (*Optional*) Specific tissue from which the cells were derived, key for context and specificity in cell studies.                                                                                                                                                                                                                                                                                                                                                                              |
| `obs["tissue_ontology_term_id"]`                  | `string`  | (*Optional*) Ontology term identifier for the tissue, providing a standardized reference for the tissue type. For organoid or tissue samples, the Uber-anatomy ontology (`UBERON:`) is used. The term ids must be a child term of `UBERON:0001062` (anatomical entity). For cell cultures, the Cell Ontology (`CL:`) is used. The term ids cannot be `CL:0000255`, `CL:0000257` or `CL:0000548`.                                                                                              |
| `obs["tissue_general"]`                           | `string`  | (*Optional*) General category or classification of the tissue, useful for broader grouping and comparison of cell data.                                                                                                                                                                                                                                                                                                                                                                       |
| `obs["tissue_general_ontology_term_id"]`          | `string`  | (*Optional*) Ontology term identifier for the general tissue category, aiding in standardizing and grouping tissue types. For organoid or tissue samples, the Uber-anatomy ontology (`UBERON:`) is used. The term ids must be a child term of `UBERON:0001062` (anatomical entity). For cell cultures, the Cell Ontology (`CL:`) is used. The term ids cannot be `CL:0000255`, `CL:0000257` or `CL:0000548`.                                                                                  |
| `obs["batch"]`                                    | `string`  | (*Optional*) A batch identifier. This label is very context-dependent and may be a combination of the tissue, assay, donor, etc.                                                                                                                                                                                                                                                                                                                                                              |
| `obs["soma_joinid"]`                              | `integer` | (*Optional*) If the dataset was retrieved from CELLxGENE census, this is a unique identifier for the cell.                                                                                                                                                                                                                                                                                                                                                                                    |
| `obs["size_factors"]`                             | `double`  | (*Optional*) The size factors created by the normalisation method, if any.                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `var["feature_id"]`                               | `string`  | (*Optional*) Unique identifier for the feature, usually a ENSEMBL gene id.                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `var["feature_name"]`                             | `string`  | A human-readable name for the feature, usually a gene symbol.                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `var["soma_joinid"]`                              | `integer` | (*Optional*) If the dataset was retrieved from CELLxGENE census, this is a unique identifier for the feature.                                                                                                                                                                                                                                                                                                                                                                                 |
| `var["hvg"]`                                      | `boolean` | Whether or not the feature is considered to be a ‘highly variable gene’.                                                                                                                                                                                                                                                                                                                                                                                                                      |
| `var["hvg_score"]`                                | `double`  | A score for the feature indicating how highly variable it is.                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `obsm["X_pca"]`                                   | `double`  | The resulting PCA embedding.                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| `obsp["knn_distances"]`                           | `double`  | K nearest neighbors distance matrix.                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| `obsp["knn_connectivities"]`                      | `double`  | K nearest neighbors connectivities matrix.                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `varm["pca_loadings"]`                            | `double`  | The PCA loadings matrix.                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| `layers["counts"]`                                | `integer` | Raw counts.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| `layers["normalized"]`                            | `double`  | Normalised expression values.                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `uns["dataset_id"]`                               | `string`  | A unique identifier for the dataset. This is different from the `obs.dataset_id` field, which is the identifier for the dataset from which the cell data is derived.                                                                                                                                                                                                                                                                                                                          |
| `uns["dataset_name"]`                             | `string`  | A human-readable name for the dataset.                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| `uns["dataset_url"]`                              | `string`  | (*Optional*) Link to the original source of the dataset.                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| `uns["dataset_reference"]`                        | `string`  | (*Optional*) Bibtex reference of the paper in which the dataset was published.                                                                                                                                                                                                                                                                                                                                                                                                                |
| `uns["dataset_summary"]`                          | `string`  | Short description of the dataset.                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| `uns["dataset_description"]`                      | `string`  | Long description of the dataset.                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| `uns["dataset_organism"]`                         | `string`  | (*Optional*) The organism of the sample in the dataset.                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| `uns["normalization_id"]`                         | `string`  | Which normalization was used.                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `uns["pca_variance"]`                             | `double`  | The PCA variance objects.                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| `uns["knn"]`                                      | `object`  | Supplementary K nearest neighbors data.                                                                                                                                                                                                                                                                                                                                                                                                                                                       |

</div>

## Component type: Data processor

Path:
[`src/dimensionality_reduction`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/dimensionality_reduction)

A dimensionality reduction dataset processor.

Arguments:

<div class="small">

| Name                | Type   | Description                                                    |
|:--------------------|:-------|:---------------------------------------------------------------|
| `--input`           | `file` | A dataset processed by the common dataset processing pipeline. |
| `--output_dataset`  | `file` | (*Output*) The dataset to pass to a method.                    |
| `--output_solution` | `file` | (*Output*) The data for evaluating a dimensionality reduction. |

</div>

## File format: Dataset

The dataset to pass to a method.

Example file:
`resources_test/dimensionality_reduction/pancreas/dataset.h5ad`

Format:

<div class="small">

    AnnData object
     var: 'hvg_score'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'normalization_id'

</div>

Slot description:

<div class="small">

| Slot                      | Type      | Description                                                                          |
|:--------------------------|:----------|:-------------------------------------------------------------------------------------|
| `var["hvg_score"]`        | `double`  | High variability gene score (normalized dispersion). The greater, the more variable. |
| `layers["counts"]`        | `integer` | Raw counts.                                                                          |
| `layers["normalized"]`    | `double`  | Normalized expression values.                                                        |
| `uns["dataset_id"]`       | `string`  | A unique identifier for the dataset.                                                 |
| `uns["normalization_id"]` | `string`  | Which normalization was used.                                                        |

</div>

## File format: Test data

The data for evaluating a dimensionality reduction.

Example file:
`resources_test/dimensionality_reduction/pancreas/solution.h5ad`

Format:

<div class="small">

    AnnData object
     obs: 'cell_type'
     var: 'hvg_score'
     layers: 'counts', 'normalized'
     uns: 'dataset_id', 'dataset_name', 'dataset_url', 'dataset_reference', 'dataset_summary', 'dataset_description', 'dataset_organism', 'normalization_id'

</div>

Slot description:

<div class="small">

| Slot                         | Type      | Description                                                                                              |
|:-----------------------------|:----------|:---------------------------------------------------------------------------------------------------------|
| `obs["cell_type"]`           | `string`  | Classification of the cell type based on its characteristics and function within the tissue or organism. |
| `var["hvg_score"]`           | `double`  | High variability gene score (normalized dispersion). The greater, the more variable.                     |
| `layers["counts"]`           | `integer` | Raw counts.                                                                                              |
| `layers["normalized"]`       | `double`  | Normalized expression values.                                                                            |
| `uns["dataset_id"]`          | `string`  | A unique identifier for the dataset.                                                                     |
| `uns["dataset_name"]`        | `string`  | Nicely formatted name.                                                                                   |
| `uns["dataset_url"]`         | `string`  | (*Optional*) Link to the original source of the dataset.                                                 |
| `uns["dataset_reference"]`   | `string`  | (*Optional*) Bibtex reference of the paper in which the dataset was published.                           |
| `uns["dataset_summary"]`     | `string`  | Short description of the dataset.                                                                        |
| `uns["dataset_description"]` | `string`  | Long description of the dataset.                                                                         |
| `uns["dataset_organism"]`    | `string`  | (*Optional*) The organism of the sample in the dataset.                                                  |
| `uns["normalization_id"]`    | `string`  | Which normalization was used.                                                                            |

</div>

## Component type: Control method

Path:
[`src/dimensionality_reduction/control_methods`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/dimensionality_reduction/control_methods)

Quality control methods for verifying the pipeline.

Arguments:

<div class="small">

| Name               | Type   | Description                                                   |
|:-------------------|:-------|:--------------------------------------------------------------|
| `--input`          | `file` | The dataset to pass to a method.                              |
| `--input_solution` | `file` | The data for evaluating a dimensionality reduction.           |
| `--output`         | `file` | (*Output*) A dataset with dimensionality reduction embedding. |

</div>

## Component type: Method

Path:
[`src/dimensionality_reduction/methods`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/dimensionality_reduction/methods)

A dimensionality reduction method.

Arguments:

<div class="small">

| Name       | Type   | Description                                                   |
|:-----------|:-------|:--------------------------------------------------------------|
| `--input`  | `file` | The dataset to pass to a method.                              |
| `--output` | `file` | (*Output*) A dataset with dimensionality reduction embedding. |

</div>

## Component type: Metric

Path:
[`src/dimensionality_reduction/metrics`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/dimensionality_reduction/metrics)

A dimensionality reduction metric.

Arguments:

<div class="small">

| Name                | Type   | Description                                         |
|:--------------------|:-------|:----------------------------------------------------|
| `--input_embedding` | `file` | A dataset with dimensionality reduction embedding.  |
| `--input_solution`  | `file` | The data for evaluating a dimensionality reduction. |
| `--output`          | `file` | (*Output*) Metric score file.                       |

</div>

## File format: Embedding

A dataset with dimensionality reduction embedding.

Example file:
`resources_test/dimensionality_reduction/pancreas/embedding.h5ad`

Format:

<div class="small">

    AnnData object
     obsm: 'X_emb'
     uns: 'dataset_id', 'method_id', 'normalization_id'

</div>

Slot description:

<div class="small">

| Slot                      | Type     | Description                          |
|:--------------------------|:---------|:-------------------------------------|
| `obsm["X_emb"]`           | `double` | The dimensionally reduced embedding. |
| `uns["dataset_id"]`       | `string` | A unique identifier for the dataset. |
| `uns["method_id"]`        | `string` | A unique identifier for the method.  |
| `uns["normalization_id"]` | `string` | Which normalization was used.        |

</div>

## File format: Score

Metric score file

Example file:
`resources_test/dimensionality_reduction/pancreas/score.h5ad`

Format:

<div class="small">

    AnnData object
     uns: 'dataset_id', 'normalization_id', 'method_id', 'metric_ids', 'metric_values'

</div>

Slot description:

<div class="small">

| Slot                      | Type     | Description                                                                                  |
|:--------------------------|:---------|:---------------------------------------------------------------------------------------------|
| `uns["dataset_id"]`       | `string` | A unique identifier for the dataset.                                                         |
| `uns["normalization_id"]` | `string` | Which normalization was used.                                                                |
| `uns["method_id"]`        | `string` | A unique identifier for the method.                                                          |
| `uns["metric_ids"]`       | `string` | One or more unique metric identifiers.                                                       |
| `uns["metric_values"]`    | `double` | The metric values obtained for the given prediction. Must be of same length as ‘metric_ids’. |

</div>

