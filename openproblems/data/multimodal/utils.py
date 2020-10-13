import numpy as np
import pandas as pd
import scprep
import anndata

from ..utils import filter_genes_cells


def filter_joint_data_empty_cells(adata):
    assert np.all(adata.uns["mode2_obs"] == adata.obs.index)
    n_mode1 = scprep.utils.toarray(adata.X.sum(axis=1)).flatten()
    n_mode2 = scprep.utils.toarray(adata.obsm["mode2"].sum(axis=1)).flatten()
    keep_cells = np.minimum(n_mode1, n_mode2) > 0
    adata.uns["mode2_obs"] = adata.uns["mode2_obs"][keep_cells]
    adata = adata[keep_cells, :]
    return adata.copy()


def create_joint_adata(
    X, Y, X_index=None, X_columns=None, Y_index=None, Y_columns=None
):
    if X_index is None:
        X_index = X.index
    if X_columns is None:
        X_columns = X.columns
    if Y_index is None:
        Y_index = Y.index
    if Y_columns is None:
        Y_columns = Y.columns
    joint_index = np.sort(np.intersect1d(X_index, Y_index))
    try:
        X = X.loc[joint_index]
        Y = Y.loc[joint_index]
    except AttributeError:
        # keep only common observations
        X_keep_idx = np.isin(X_index, joint_index)
        Y_keep_idx = np.isin(Y_index, joint_index)
        X = X[X_keep_idx]
        Y = Y[Y_keep_idx]

        # reorder by alphabetical
        X_index_sub = scprep.utils.toarray(X_index[X_keep_idx])
        Y_index_sub = scprep.utils.toarray(Y_index[Y_keep_idx])
        X = X[np.argsort(X_index_sub)]
        Y = Y[np.argsort(Y_index_sub)]

        # check order is correct
        assert (X_index_sub[np.argsort(X_index_sub)] == joint_index).all()
        assert (Y_index_sub[np.argsort(Y_index_sub)] == joint_index).all()
    adata = anndata.AnnData(
        scprep.utils.to_array_or_spmatrix(X).tocsr(),
        obs=pd.DataFrame(index=joint_index),
        var=pd.DataFrame(index=X_columns),
    )
    adata.obsm["mode2"] = scprep.utils.to_array_or_spmatrix(Y).tocsr()
    adata.uns["mode2_obs"] = joint_index
    adata.uns["mode2_var"] = scprep.utils.toarray(Y_columns)
    return adata


def subset_joint_data(adata, n_cells=500, n_genes=200):
    keep_cells = np.random.choice(adata.shape[0], n_cells, replace=False)
    adata = adata[keep_cells]
    adata.uns["mode2_obs"] = adata.uns["mode2_obs"][keep_cells]

    keep_rna_genes = np.argwhere(
        scprep.utils.toarray(adata.X.sum(axis=0)).flatten() > 0
    ).flatten()

    if len(keep_rna_genes) > n_genes:
        keep_rna_genes = keep_rna_genes[
            np.random.choice(len(keep_rna_genes), n_genes, replace=False)
        ]
    adata = adata[:, keep_rna_genes].copy()

    keep_adt_genes = np.argwhere(
        scprep.utils.toarray(adata.obsm["mode2"].sum(axis=0)).flatten() > 0
    ).flatten()
    if len(keep_adt_genes) > n_genes:
        keep_adt_genes = keep_adt_genes[
            np.random.choice(len(keep_adt_genes), n_genes, replace=False)
        ]
    adata.obsm["mode2"] = adata.obsm["mode2"][:, keep_adt_genes]
    adata.uns["mode2_var"] = adata.uns["mode2_var"][keep_adt_genes]

    adata = filter_joint_data_empty_cells(adata)
    return adata
