from anndata import AnnData

import numpy as np


def obs_means(adata: AnnData, cluster_key: str) -> AnnData:
    """Return means over observation key."""

    labels = adata.obs[cluster_key].cat.categories
    means = np.empty((labels.shape[0], adata.shape[1]))
    for i, lab in enumerate(labels):
        means[i, :] = adata[adata.obs[cluster_key] == lab].X.mean(axis=0).flatten()
    adata_means = AnnData(means)
    adata_means.obs_names = labels
    adata_means.var_names = adata.var_names

    return adata_means


def normalize_coefficients(_prop: np.array) -> np.array:
    """Normalize coefficients to sum to 1."""
    prop = _prop.copy()
    prop[prop < 0] = 0
    prop = prop / prop.sum(axis=1, keepdims=1)
    return prop
