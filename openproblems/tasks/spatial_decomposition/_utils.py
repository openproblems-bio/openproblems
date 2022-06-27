from typing import Tuple

import anndata as ad
import numpy as np
import pandas as pd


def merge_sc_and_sp(
    adata_sc: ad.AnnData,
    adata_sp: ad.AnnData,
    batch_key: str = "modality",
) -> ad.AnnData:

    for k in adata_sp.obsm.keys():
        if isinstance(adata_sp.obsm[k], pd.DataFrame):
            n_col = adata_sp.obsm[k].shape[1]
            n_row = adata_sc.shape[0]
            df = pd.DataFrame(
                np.zeros((n_row, n_col)),
                columns=adata_sp.obsm[k].columns.copy(),
                index=adata_sc.obs_names.copy(),
            )
            adata_sc.obsm[k] = df

    # merge single cell and spatial data
    adata_merged = ad.concat(
        {"sp": adata_sp, "sc": adata_sc},
        label=batch_key,
        join="outer",
        index_unique=None,
        merge="unique",
        uns_merge="unique",
    )

    adata_merged.obs["label"] = pd.Categorical(adata_merged.obs["label"])

    return adata_merged


def split_sc_and_sp(
    adata_merged: ad.AnnData,
    batch_key: str = "modality",
) -> Tuple[ad.AnnData, ad.AnnData]:

    # split single cell and spatial data
    is_sp = adata_merged.obs[batch_key] == "sp"
    adata_sp = adata_merged[is_sp, :].copy()
    adata_sc = adata_merged[~is_sp, :].copy()

    return (adata_sc, adata_sp)


def obs_means(adata: ad.AnnData, cluster_key: str) -> ad.AnnData:
    """Return means over observation key."""

    labels = adata.obs[cluster_key].cat.categories
    means = np.empty((labels.shape[0], adata.shape[1]))
    for i, lab in enumerate(labels):
        means[i, :] = adata[adata.obs[cluster_key] == lab].X.mean(axis=0).flatten()
    adata_means = ad.AnnData(means)
    adata_means.obs_names = labels
    adata_means.var_names = adata.var_names

    return adata_means


def normalize_coefficients(_prop: np.array) -> np.array:
    """Normalize coefficients to sum to 1."""
    prop = _prop.copy()
    prop[prop < 0] = 0
    prop = prop / prop.sum(axis=1, keepdims=1)
    return prop
