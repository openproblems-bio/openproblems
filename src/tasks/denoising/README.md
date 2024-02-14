# Denoising

Removing noise in sparse single-cell RNA-sequencing count data

Path:
[`src/tasks/denoising`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/tasks/denoising)

## Motivation

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

## Description

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

## Authors & contributors

| name              | roles              |
|:------------------|:-------------------|
| Wesley Lewis      | author, maintainer |
| Scott Gigante     | author, maintainer |
| Robrecht Cannoodt | author             |
| Kai Waldrant      | author             |

## API

``` mermaid
flowchart LR
  file_common_dataset("Common dataset")
  comp_process_dataset[/"Data processor"/]
  file_train("Training data")
  file_test("Test data")
  comp_control_method[/"Control method"/]
  comp_method[/"Method"/]
  comp_metric[/"Metric"/]
  file_denoised("Denoised data")
  file_score("Score")
  file_common_dataset---comp_process_dataset
  comp_process_dataset-->file_train
  comp_process_dataset-->file_test
  file_train---comp_control_method
  file_train---comp_method
  file_test---comp_control_method
  file_test---comp_metric
  comp_control_method-->file_denoised
  comp_method-->file_denoised
  comp_metric-->file_score
  file_denoised---comp_metric
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
| `var["hvg_score"]`                                | `integer` | A ranking of the features by hvg.                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
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
[`src/denoising`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/denoising)

A denoising dataset processor.

Arguments:

<div class="small">

| Name             | Type   | Description                                                       |
|:-----------------|:-------|:------------------------------------------------------------------|
| `--input`        | `file` | A dataset processed by the common dataset processing pipeline.    |
| `--output_train` | `file` | (*Output*) The subset of molecules used for the training dataset. |
| `--output_test`  | `file` | (*Output*) The subset of molecules used for the test dataset.     |

</div>

## File format: Training data

The subset of molecules used for the training dataset

Example file: `resources_test/denoising/pancreas/train.h5ad`

Format:

<div class="small">

    AnnData object
     layers: 'counts'
     uns: 'dataset_id'

</div>

Slot description:

<div class="small">

| Slot                | Type      | Description                          |
|:--------------------|:----------|:-------------------------------------|
| `layers["counts"]`  | `integer` | Raw counts.                          |
| `uns["dataset_id"]` | `string`  | A unique identifier for the dataset. |

</div>

## File format: Test data

The subset of molecules used for the test dataset

Example file: `resources_test/denoising/pancreas/test.h5ad`

Format:

<div class="small">

    AnnData object
     layers: 'counts'
     uns: 'dataset_id', 'dataset_name', 'dataset_url', 'dataset_reference', 'dataset_summary', 'dataset_description', 'dataset_organism'

</div>

Slot description:

<div class="small">

| Slot                         | Type      | Description                                                                    |
|:-----------------------------|:----------|:-------------------------------------------------------------------------------|
| `layers["counts"]`           | `integer` | Raw counts.                                                                    |
| `uns["dataset_id"]`          | `string`  | A unique identifier for the dataset.                                           |
| `uns["dataset_name"]`        | `string`  | Nicely formatted name.                                                         |
| `uns["dataset_url"]`         | `string`  | (*Optional*) Link to the original source of the dataset.                       |
| `uns["dataset_reference"]`   | `string`  | (*Optional*) Bibtex reference of the paper in which the dataset was published. |
| `uns["dataset_summary"]`     | `string`  | Short description of the dataset.                                              |
| `uns["dataset_description"]` | `string`  | Long description of the dataset.                                               |
| `uns["dataset_organism"]`    | `string`  | (*Optional*) The organism of the sample in the dataset.                        |

</div>

## Component type: Control method

Path:
[`src/denoising/control_methods`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/denoising/control_methods)

Quality control methods for verifying the pipeline.

Arguments:

<div class="small">

| Name            | Type   | Description                                                    |
|:----------------|:-------|:---------------------------------------------------------------|
| `--input_train` | `file` | The subset of molecules used for the training dataset.         |
| `--input_test`  | `file` | The subset of molecules used for the test dataset.             |
| `--output`      | `file` | (*Output*) A denoised dataset as output by a denoising method. |

</div>

## Component type: Method

Path:
[`src/denoising/methods`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/denoising/methods)

A denoising method.

Arguments:

<div class="small">

| Name            | Type   | Description                                                    |
|:----------------|:-------|:---------------------------------------------------------------|
| `--input_train` | `file` | The subset of molecules used for the training dataset.         |
| `--output`      | `file` | (*Output*) A denoised dataset as output by a denoising method. |

</div>

## Component type: Metric

Path:
[`src/denoising/metrics`](https://github.com/openproblems-bio/openproblems-v2/tree/main/src/denoising/metrics)

A denoising metric.

Arguments:

<div class="small">

| Name               | Type   | Description                                         |
|:-------------------|:-------|:----------------------------------------------------|
| `--input_test`     | `file` | The subset of molecules used for the test dataset.  |
| `--input_denoised` | `file` | A denoised dataset as output by a denoising method. |
| `--output`         | `file` | (*Output*) Metric score file.                       |

</div>

## File format: Denoised data

A denoised dataset as output by a denoising method.

Example file: `resources_test/denoising/pancreas/denoised.h5ad`

Format:

<div class="small">

    AnnData object
     layers: 'counts', 'denoised'
     uns: 'dataset_id', 'method_id'

</div>

Slot description:

<div class="small">

| Slot                 | Type      | Description                          |
|:---------------------|:----------|:-------------------------------------|
| `layers["counts"]`   | `integer` | Raw counts.                          |
| `layers["denoised"]` | `integer` | denoised data.                       |
| `uns["dataset_id"]`  | `string`  | A unique identifier for the dataset. |
| `uns["method_id"]`   | `string`  | A unique identifier for the method.  |

</div>

## File format: Score

NA

Example file: `resources_test/denoising/pancreas/score.h5ad`

Description:

Metric score file

Format:

<div class="small">

    AnnData object
     uns: 'dataset_id', 'method_id', 'metric_ids', 'metric_values'

</div>

Slot description:

<div class="small">

| Slot                   | Type     | Description                                                                                  |
|:-----------------------|:---------|:---------------------------------------------------------------------------------------------|
| `uns["dataset_id"]`    | `string` | A unique identifier for the dataset.                                                         |
| `uns["method_id"]`     | `string` | A unique identifier for the method.                                                          |
| `uns["metric_ids"]`    | `string` | One or more unique metric identifiers.                                                       |
| `uns["metric_values"]` | `double` | The metric values obtained for the given prediction. Must be of same length as ‘metric_ids’. |

</div>
