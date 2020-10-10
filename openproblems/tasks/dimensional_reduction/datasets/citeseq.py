import os
import tempfile

import numpy as np
import pandas as pd
import scanpy as sc
import scprep
import anndata

from .utils import loader

ADT_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format=file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DADT%5Fumi%2Ecsv%2Egz"
RNA_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format=file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DRNA%5Fumi%2Ecsv%2Egz"


@loader
def load_citeseq_cbmc(test=False):
    if test:
        adata = load_citeseq_cbmc(test=False)
        keep_cells = np.random.choice(adata.shape[0], 500, replace=False)
        adata = adata[keep_cells]

        keep_rna_genes = np.argwhere(
            scprep.utils.toarray(adata.X.sum(axis=0)).flatten() > 0
        ).flatten()

        if len(keep_rna_genes) > 200:
            keep_rna_genes = keep_rna_genes[
                np.random.choice(len(keep_rna_genes), 200, replace=False)
            ]
        adata = adata[:, keep_rna_genes].copy()
        adata.uns["mode2_obs"] = adata.uns["mode2_obs"][keep_cells]

        keep_adt_genes = np.argwhere(
            scprep.utils.toarray(adata.obsm["mode2"].sum(axis=0)).flatten() > 0
        ).flatten()
        if len(keep_adt_genes) > 200:
            keep_adt_genes = keep_adt_genes[
                np.random.choice(len(keep_adt_genes), 200, replace=False)
            ]
        adata.obsm["mode2"] = adata.obsm["mode2"][:, keep_adt_genes]
        adata.uns["mode2_var"] = adata.uns["mode2_var"][keep_adt_genes]
        return adata

    rna_data = scprep.io.load_csv(
        RNA_URL, cell_axis="col", compression="gzip", sparse=True
    )
    adt_data = scprep.io.load_csv(
        ADT_URL, cell_axis="col", compression="gzip", sparse=True
    )

    rna_data = scprep.filter.filter_empty_cells(rna_data)
    adt_data = scprep.filter.filter_empty_cells(adt_data)

    common_cells = np.intersect1d(rna_data.index, adt_data.index)

    rna_data = rna_data.loc[common_cells]
    adt_data = adt_data.loc[common_cells]

    adata = anndata.AnnData(
        rna_data.sparse.to_coo().tocsr(),
        obs=pd.DataFrame(index=rna_data.index),
        var=pd.DataFrame(index=rna_data.columns),
    )
    adata.obsm["mode2"] = adt_data.sparse.to_coo().tocsr()
    adata.uns["mode2_obs"] = adt_data.index.to_numpy()
    adata.uns["mode2_var"] = adt_data.columns.to_numpy()

    return adata
