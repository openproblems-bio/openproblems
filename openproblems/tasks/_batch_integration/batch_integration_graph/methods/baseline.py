from .....tools.decorators import baseline_method
from .....tools.utils import check_version

import numpy as np
import scanpy as sc


def _set_uns(adata):
    adata.uns["neighbors"] = adata.uns["uni"]
    adata.uns["neighbors"]["connectivities_key"] = "connectivities"
    adata.uns["neighbors"]["distances_key"] = "distances"


@baseline_method(method_name="No Integration", method_summary="TODO")
def no_integration(adata, test=False):
    adata.obsp["connectivities"] = adata.obsp["uni_connectivities"]
    adata.obsp["distances"] = adata.obsp["uni_distances"]
    _set_uns(adata)
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


def _randomize_features(X, partition=None):
    X_out = X.copy()
    if partition is None:
        partition = np.full(X.shape[0], 0)
    else:
        partition = np.asarray(partition)
    for partition_name in np.unique(partition):
        partition_idx = np.argwhere(partition == partition_name).flatten()
        X_out[partition_idx] = X[np.random.permutation(partition_idx)]
    return X_out


def _randomize_graph(adata, partition=None):
    distances, connectivities = (
        adata.obsp["uni_distances"],
        adata.obsp["uni_connectivities"],
    )
    new_idx = _randomize_features(np.arange(distances.shape[0]), partition=partition)
    adata.obsp["distances"] = distances[new_idx][:, new_idx]
    adata.obsp["connectivities"] = connectivities[new_idx][:, new_idx]
    _set_uns(adata)
    return adata


def _random_embedding(partition):
    from sklearn.preprocessing import LabelEncoder
    from sklearn.preprocessing import OneHotEncoder

    embedding = OneHotEncoder().fit_transform(
        LabelEncoder().fit_transform(partition)[:, None]
    )
    embedding = embedding + np.random.uniform(-0.01, 0.01, embedding.shape)
    return embedding


@baseline_method(method_name="Random Integration", method_summary="TODO")
def random_integration(adata, test=False):
    adata = _randomize_graph(adata)
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(method_name="Random Integration by Celltype", method_summary="TODO")
def celltype_random_integration(adata, test=False):
    adata = _randomize_graph(
        adata,
        partition=adata.obs["labels"].to_numpy(),
    )
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(method_name="Random Integration by Batch", method_summary="TODO")
def batch_random_integration(adata, test=False):
    adata = _randomize_graph(
        adata,
        partition=adata.obs["batch"].to_numpy(),
    )
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(method_name="Random Graph by Celltype", method_summary="TODO")
def celltype_random_graph(adata, test=False):
    adata.obsm["X_emb"] = _random_embedding(partition=adata.obs["labels"])
    sc.pp.neighbors(adata, use_rep="X_emb")
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata
