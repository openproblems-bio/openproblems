# Spatial Decomposition/Deconvolution

## The task

Spatial decomposition (also often referred to as Spatial deconvolution) is
applicable to spatial transcriptomics data where the transcription profile of
each capture location (spot, voxel, bead, etc.) do not share a bijective
relationship with the cells in the tissue, i.e., multiple cells may contribute
to the same capture location. The task of spatial decomposition then refers to
estimating the composition of cell types/states that are present at each capture
location. The cell type/states estimates are presented as proportion values,
representing the proportion of the cells at each capture location that belong to
a given cell type.

We distinguish between _reference-based_ decomposition and _de novo_
decomposition, where the former leverage external data (e.g., scRNA-seq or
scNuc-seq) to guide the inference process, while the latter only work with the
spatial data. We require that all datasets have an associated reference single
cell data set, but methods are free to ignore this information.

## Metrics

### R2

R2 pronounced as "R squared", also known as the "coefficient of determination". R2
reports the fraction of the true proportion values' (`adata.obsm["proportions_true"]`)
variance that can be explained by the predicted proportion values
(`adata.obsm["proportion_pred"]`). The **best score**, and upper bound, is 1.0. There is
no fixed lower bound for the metric. The _uniform/non-weighted average_ across all cell
types/states is used to summarize performance. See the
[sklearn](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.r2_score.html)
documentation for details on the implementation and the
[wikipedia](https://en.wikipedia.org/wiki/Coefficient_of_determination) site for more
general information regarding the metric.

## API

Datasets consists of 2 `anndata.AnnData` objects, concatenated by key
`adata.obs["modality"]` with values:

* `sc` for the single cell reference.
* `sp` for the target spatial dataset.

In the single cell reference, cell-types are stored in `adata_sc.obs["label"]`. No label
may have fewer than 25 cells.

In the spatial target, ground-truth cell-type proportions are stored in
`adata_spatial.obsm["proportions_true"]`.

Methods should return only the spatial data with inferred proportions stored in
`adata_spatial.obsm["proportions_pred"]`.

Metrics shall compare `adata_spatial.obsm['proportions_pred']` to
`adata_spatial.obsm['proportions_true']`.
