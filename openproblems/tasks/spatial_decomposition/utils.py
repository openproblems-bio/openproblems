from typing import Optional
from typing import Tuple

import anndata as ad
import numpy as np


def merge_sc_and_sp(
    adata_sc: ad.AnnData,
    adata_sp: ad.AnnData,
    batch_key: str = "modality",
    test: bool = False,
    test_n_genes: int = 1000,
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

    if test:
        adata_merged = adata_merged[:, :test_n_genes].copy()

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


def obs_means(
    adata: ad.AnnData, cluster_key: str, obsm: Optional[str] = None
) -> ad.AnnData:
    """Return means over observation key."""

    labels = adata.obs[cluster_key].cat.categories
    n_var = adata.shape[1] if obsm is None else adata.obsm[obsm].shape[1]
    means = np.empty((labels.shape[0], n_var))
    for i, lab in enumerate(labels):
        adata_lab = adata[adata.obs[cluster_key] == lab]
        x_lab = adata_lab.X if obsm is None else adata_lab.obsm[obsm]
        means[i, :] = x_lab.mean(axis=0).flatten()
    adata_means = ad.AnnData(means)
    adata_means.obs_names = labels
    if obsm is None:
        adata_means.var_names = adata.var_names

    return adata_means


def normalize_coefficients(prop: np.array) -> np.array:
    """Normalize coefficients to sum to 1."""
    prop = prop.copy()
    prop[prop < 0] = 0
    prop = prop / prop.sum(axis=1, keepdims=1)
    return prop
