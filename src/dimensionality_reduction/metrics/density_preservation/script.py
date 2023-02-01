

import anndata as ad
import numpy as np
from typing import Optional
from umap import UMAP
from scipy.stats import pearsonr

## VIASH START
par = {
    "input_reduced": "resources_test/dimensionality_reduction/pancreas/reduced.h5ad",
    "input_test": "resources_test/dimensionality_reduction/pancreas/test.h5ad",
    "output": "score.h5ad",
}
## VIASH END

# Interpreted from:
# https://github.com/lmcinnes/umap/blob/317ce81dc64aec9e279aa1374ac809d9ced236f6/umap/umap_.py#L1190-L1243
#
# Author: Leland McInnes <leland.mcinnes@gmail.com>
#
# License: BSD 3 clause
def _calculate_radii(
    X: np.ndarray,
    n_neighbors: int = 30,
    random_state: Optional[int] = None
) -> np.ndarray:
    from umap.umap_ import fuzzy_simplicial_set
    from umap.umap_ import nearest_neighbors

    (knn_indices, knn_dists, _) = nearest_neighbors(
        X,
        n_neighbors,
        "euclidean",
        {},
        False,
        random_state,
        verbose=False,
    )

    emb_graph, _, _, emb_dists = fuzzy_simplicial_set(
        X,
        n_neighbors,
        random_state,
        "euclidean",
        {},
        knn_indices,
        knn_dists,
        verbose=False,
        return_dists=True,
    )

    emb_graph = emb_graph.tocoo()
    emb_graph.sum_duplicates()
    emb_graph.eliminate_zeros()

    n_vertices = emb_graph.shape[1]

    mu_sum = np.zeros(n_vertices, dtype=np.float32)
    re = np.zeros(n_vertices, dtype=np.float32)

    head = emb_graph.row
    tail = emb_graph.col
    for i in range(len(head)):
        j = head[i]
        k = tail[i]
        D = emb_dists[j, k]
        mu = emb_graph.data[i]
        re[j] += mu * D
        re[k] += mu * D
        mu_sum[j] += mu
        mu_sum[k] += mu

    epsilon = 1e-8
    return np.log(epsilon + (re / mu_sum))

def compute_density_preservation(
    X_emb: np.ndarray,
    high_dim: np.ndarray,
    n_neighbors: int = 30,
    random_state: Optional[int] = None
) -> float:
    if np.any(np.isnan(X_emb)):
        return 0.0
    
    print("Compute local radii in original data", flush=True)
    _, ro, _ = UMAP(
        n_neighbors=_K,
        random_state=_SEED,
        densmap=True,
        output_dens=True
    ).fit_transform(high_dim)

    print("Compute local radii of embedding", flush=True)
    re = _calculate_radii(
        X_emb,
        n_neighbors=_K,
        random_state=_SEED
    )
    
    print("Compute pearson correlation", flush=True)
    return pearsonr(ro, re)[0]

# number of neighbors
_K = 30
# Fix seed
_SEED = 42

print("Load data", flush=True)
input_test = ad.read_h5ad(par["input_test"])
input_reduced = ad.read_h5ad(par["input_reduced"])

high_dim = input_test.layers["normalized"]
X_emb = input_reduced.obsm["X_emb"]

density_preservation = compute_density_preservation(
    X_emb=X_emb,
    high_dim=high_dim,
    n_neighbors=_K,
    random_state=_SEED
)

print("Create output AnnData object", flush=True)
output = ad.AnnData(
    uns={
        "dataset_id": input_test.uns["dataset_id"],
        "normalization_id": input_test.uns["normalization_id"],
        "method_id": input_reduced.uns["method_id"],
        "metric_ids": [ "density_preservation" ],
        "metric_values": [ density_preservation ]
    }
)

print("Write data to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")