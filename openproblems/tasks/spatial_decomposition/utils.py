from typing import Tuple

import anndata as ad
import numpy as np


def merge_sc_and_sp(
    adata_sc: ad.AnnData,
    adata_sp: ad.AnnData,
    batch_key: str = "modality",
) -> ad.AnnData:

    # merge single cell and spatial data
    adata_merged = ad.concat(
        {"sp": adata_sp, "sc": adata_sc},
        label=batch_key,
        join="outer",
        index_unique=None,
        merge="unique",
        uns_merge="unique",
    )

    adata_merged.strings_to_categoricals()

    return adata_merged


def split_sc_and_sp(
    adata_merged: ad.AnnData,
    batch_key: str = "modality",
) -> Tuple[ad.AnnData, ad.AnnData]:

    # split single cell and spatial data
    is_sp = adata_merged.obs[batch_key] == "sp"
    adata_sp = adata_merged[is_sp, :].copy()
    adata_sc = adata_merged[~is_sp, :].copy()

    return adata_sc, adata_sp


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


def normalize_coefficients(prop: np.array) -> np.array:
    """Normalize coefficients to sum to 1."""
    prop = prop.copy()
    prop[prop < 0] = 0
    prop = prop / prop.sum(axis=1, keepdims=1)
    return prop
