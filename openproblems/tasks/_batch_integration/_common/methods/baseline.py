from .....tools.decorators import baseline_method
from .....tools.utils import check_version

import numpy as np


def _set_uns(adata):
    adata.uns["neighbors"] = adata.uns["uni"]
    adata.uns["neighbors"]["connectivities_key"] = "connectivities"
    adata.uns["neighbors"]["distances_key"] = "distances"


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


def _random_embedding(partition, jitter=0.01):
    from sklearn.preprocessing import LabelEncoder
    from sklearn.preprocessing import OneHotEncoder

    embedding = OneHotEncoder().fit_transform(
        LabelEncoder().fit_transform(partition)[:, None]
    )
    if jitter is not None:
        embedding = embedding + np.random.uniform(-1 * jitter, jitter, embedding.shape)
    return embedding


@baseline_method(
    method_name="No Integration",
    method_summary=(
        "Cells are embedded by PCA on the unintegrated data. A graph is built on this"
        " PCA embedding."
    ),
)
def no_integration(adata, test=False):
    adata.obsp["connectivities"] = adata.obsp["uni_connectivities"]
    adata.obsp["distances"] = adata.obsp["uni_distances"]
    _set_uns(adata)
    adata.obsm["X_emb"] = adata.obsm["X_uni_pca"]
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(
    method_name="Random Integration",
    method_summary=(
        "Feature values, embedding coordinates, and graph connectivity are all randomly"
        " permuted"
    ),
)
def random_integration(adata, test=False):
    adata.X = _randomize_features(adata.X)
    adata.obsm["X_emb"] = _randomize_features(adata.obsm["X_uni_pca"])
    adata = _randomize_graph(adata)
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(
    method_name="Random Integration by Celltype",
    method_summary=(
        "Feature values, embedding coordinates, and graph connectivity are all randomly"
        " permuted within each celltype label"
    ),
)
def celltype_random_integration(adata, test=False):
    adata.obsm["X_emb"] = _randomize_features(
        adata.obsm["X_uni_pca"], partition=adata.obs["labels"]
    )
    adata.X = _randomize_features(adata.X, partition=adata.obs["labels"])
    adata = _randomize_graph(
        adata,
        partition=adata.obs["labels"].to_numpy(),
    )
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata


@baseline_method(
    method_name="Random Integration by Batch",
    method_summary=(
        "Feature values, embedding coordinates, and graph connectivity are all randomly"
        " permuted within each batch label"
    ),
)
def batch_random_integration(adata, test=False):
    adata.obsm["X_emb"] = _randomize_features(
        adata.obsm["X_uni_pca"], partition=adata.obs["batch"]
    )
    adata.X = _randomize_features(adata.X, partition=adata.obs["batch"])
    adata = _randomize_graph(
        adata,
        partition=adata.obs["batch"].to_numpy(),
    )
    adata.uns["method_code_version"] = check_version("openproblems")
    return adata
