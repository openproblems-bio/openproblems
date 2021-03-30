# Spatial Decomposition/Deconvolution

Spatial decomposition (also often referred to as Spatial deconvolution) is
applicable to spatial transcriptomics data where the transcription profile of
each capture location (spot, voxel, bead, etc.) do not share a bijective
relationship with the cells in the tissue, i.e., multiple cells may contribute
to the same capture location. The task of spatial decomposition then refers to
estimating the composition of cell types/states that are present at each capture
location. The cell type/states estimates are presented as proportion values,
representing the proportion of the cells at each capture location that belong to
a given cell type.


We distinguish between _reference-based_ decomposition and _uninformed_
decomposition, where the former leverage external data (e.g., scRNA-seq or
scNuc-seq) to guide the inference process, while the latter only work with the
spatial data.


## API

Datasets should contain the following attributes:

### Reference-based


* `adata.obsm["proportions_true"]` ground truth proportion estimates 
* `adata.obsm["spatial"]` array with spatial coordinates (x,y)
* `adata.uns["sc_reference"]` an `anndata` object that holds the single cell reference data
* `adata.uns["sc_reference"].obs["label"]` cell type/state/cluster of associated cell


Note that the entries of the label column in the single cell `anndata` object must agree with the columns of the columns in the `obsm["proportions_true"]` attribute of the spatial data.

Methods should store estimates of the cell type proportions in `adata.obsm["proportions_pred"]` 

Single cell annotations are found in `adata.uns["sc_reference"].obs['labels']`.

Metrics shall compare `adata.obsm['proportions_pred']` to `adata.obsm['proportions_true']`.
