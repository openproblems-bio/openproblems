import anndata
import numpy as np
import pandas as pd
import scprep


def subset_mode2_genes(adata, keep_genes):
    """Randomly subset genes from adata.obsm["mode2"]."""
    adata.obsm["mode2"] = adata.obsm["mode2"][:, keep_genes]
    adata.uns["mode2_var"] = adata.uns["mode2_var"][keep_genes]
    if "mode2_varnames" in adata.uns:
        for varname in adata.uns["mode2_varnames"]:
            adata.uns[varname] = adata.uns[varname][keep_genes]
    return adata


def filter_joint_data_empty_cells(adata):
    """Remove empty cells and genes from a multimodal dataset."""
    import scanpy as sc

    assert np.all(adata.uns["mode2_obs"] == adata.obs.index)
    # filter cells
    n_cells_mode1 = scprep.utils.toarray(adata.X.sum(axis=1)).flatten()
    n_cells_mode2 = scprep.utils.toarray(adata.obsm["mode2"].sum(axis=1)).flatten()
    keep_cells = np.minimum(n_cells_mode1, n_cells_mode2) > 1
    adata.uns["mode2_obs"] = adata.uns["mode2_obs"][keep_cells]
    adata = adata[keep_cells, :].copy()
    # filter genes
    sc.pp.filter_genes(adata, min_counts=1)
    n_genes_mode2 = scprep.utils.toarray(adata.obsm["mode2"].sum(axis=0)).flatten()
    keep_genes_mode2 = n_genes_mode2 > 0
    adata = subset_mode2_genes(adata, keep_genes_mode2)
    return adata


def create_joint_adata(
    X, Y, X_index=None, X_columns=None, Y_index=None, Y_columns=None
):
    """Create a multimodal dataset."""
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


def subset_joint_data(adata, n_cells=600, n_genes=1500):
    """Randomly subset a multimodal dataset."""
    if adata.shape[0] > n_cells:
        keep_cells = np.random.choice(adata.shape[0], n_cells, replace=False)
        adata = adata[keep_cells].copy()
        adata.uns["mode2_obs"] = adata.uns["mode2_obs"][keep_cells]
        adata = filter_joint_data_empty_cells(adata)

    if adata.shape[1] > n_genes:
        keep_mode1_genes = np.random.choice(adata.shape[1], n_genes, replace=False)
        adata = adata[:, keep_mode1_genes].copy()

    if adata.obsm["mode2"].shape[1] > n_genes:
        keep_genes_mode2 = np.random.choice(
            adata.obsm["mode2"].shape[1], n_genes, replace=False
        )
        adata = subset_mode2_genes(adata, keep_genes_mode2)

    adata = filter_joint_data_empty_cells(adata)
    return adata
