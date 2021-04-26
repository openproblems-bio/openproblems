# Modality alignment

Modality alignment refers to the task of combining together two datasets of different modalities of measurements (e.g., single-cell RNA sequencing and single-cell ATAC sequencing) on different observations of the same biological system. Integrating such measurements allows us to analyze the interaction between the different modalities, without requiring an explicitly joint measurement like [sci-CAR](https://doi.org/10.1126/science.aau0730) or [CITE-seq](https://doi.org/10.1038/nmeth.4380).

## API

Datasets should include matched measurements from two modalities, which are contained in `adata` and `adata.obsm["mode2"]`. The task is to align these two modalities as closely as possible, without using the known bijection between the datasets. The dataset identifier should be stored in `adata.uns["dataset_id"]`.

Methods should create joint matrices `adata.obsm["aligned"]` and `adata.obsm["mode2_aligned"]` which reside in a joint space. The method identifier should be stored in `adata.uns["method_id"]`.

Metrics should evaluate how well the cells which are known to be equivalent are aligned in the joint space. The metric identifier should be stored in `adata.uns["metric_id"]`. The metric value should be stored in `adata.uns["metric_value"]`.
