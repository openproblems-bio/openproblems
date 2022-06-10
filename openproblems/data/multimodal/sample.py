from ..utils import loader
from .scicar.mouse_kidney import atac_cells_url
from .scicar.mouse_kidney import atac_genes_url
from .scicar.mouse_kidney import rna_cells_url
from .scicar.mouse_kidney import rna_genes_url
from .utils import create_joint_adata
from .utils import filter_joint_data_empty_cells

import numpy as np
import pandas as pd
import scipy.sparse


@loader(data_url="https://openproblems.bio")
def load_sample_data(test=True):
    """Create a simple dataset to use for testing in multimodal applications."""
    assert test

    rna_genes = pd.read_csv(rna_genes_url, low_memory=False, index_col=0)
    atac_genes = pd.read_csv(atac_genes_url, low_memory=False, index_col=1)
    rna_cells = pd.read_csv(rna_cells_url, low_memory=False, index_col=0)
    atac_cells = pd.read_csv(atac_cells_url, low_memory=False, index_col=0)

    keep_cells = np.intersect1d(rna_cells.index, atac_cells.index)[:200]
    rna_cells = rna_cells.loc[keep_cells]
    atac_cells = atac_cells.loc[keep_cells]

    rna_data = scipy.sparse.csr_matrix((len(keep_cells), len(rna_genes)))
    atac_data = scipy.sparse.csr_matrix((len(keep_cells), len(atac_genes)))

    adata = create_joint_adata(
        rna_data,
        atac_data,
        X_index=rna_cells.index,
        X_columns=rna_genes.index,
        Y_index=atac_cells.index,
        Y_columns=atac_genes.index,
    )

    # merge obs and var
    adata.obs = rna_cells.loc[adata.obs.index]
    adata.var = rna_genes
    for key in atac_cells.columns:
        adata.obs[key] = atac_cells[key]
    adata.uns["mode2_varnames"] = []
    for key in atac_genes.columns:
        varname = "mode2_var_{}".format(key)
        adata.uns[varname] = atac_genes[key].values
        adata.uns["mode2_varnames"].append(varname)

    # adata = subset_joint_data(adata, n_cells=200, n_genes=5000)

    adata.X = scipy.sparse.csr_matrix(np.random.poisson(0.1, adata.X.shape)).astype(
        np.float64
    )
    adata.obsm["mode2"] = scipy.sparse.csr_matrix(
        np.random.poisson(0.1, adata.obsm["mode2"].shape)
    ).astype(np.float64)

    adata = filter_joint_data_empty_cells(adata)
    return adata
