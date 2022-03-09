# ATAC signal denoising

Chromatin accessibility data obtained from Assay for Transposase-Accessible
Chromatin using sequencing (ATAC) is very space in single cells.
To enable more robust computations (e.g. expression change predictions)
on the ATAC signal a stabalization, denoising or aggregation is required.

This is a subtask of `regulatory_effect_prediction` and a lot is copied from
there.

## API

Datasets should contain the following attributes:

* `adata.uns['species']` (ensembl species name, e.g. `"mus_musculus"`)
* `adata.uns['release']` (ensembl release, e.g. `"100"`)
* `adata.uns['mode2_var_chr']` (single cell ATAC-seq peak chromosome)
* `adata.uns['mode2_var_start']` (single cell ATAC-seq peak start position)
* `adata.uns['mode2_var_end']` (single cell ATAC-seq peak end position)
* `adata.obsm['mode2']` (cell by peak matrix of single cell ATAC-seq, which is
the ground truth)
* `adata.obsm['mode2_noisy']` (`mode2` with more noise and same observations
and features)
* `adata.X` (cell by gene matrix of single cell gene expression)

## Method output
Methods shoud write denoise `adata.obsm['mode2']` and write the result into
`adata.obsm['mode2_denoised']`.
Metrics function should compare `adata.obsm['mode2_denoised']` with the true
ATAC signal in `adata.obsm['mode2']`, and return float metrics.
