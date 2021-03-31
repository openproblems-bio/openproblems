from anndata import AnnData
from openproblems.tools.decorators import metric
from scipy.sparse import issparse
from scipy.stats import pearsonr
from typing import Optional

import numpy as np


def _calculate_radii(
    X: np.ndarray, n_neighbors: int = 30, random_state: Optional[int] = None
) -> np.ndarray:
    from umap.umap_ import fuzzy_simplicial_set
    from umap.umap_ import nearest_neighbors

    # directly taken from: https://github.com/lmcinnes/umap/blob/
    # 317ce81dc64aec9e279aa1374ac809d9ced236f6/umap/umap_.py#L1190-L1243
    (knn_indices, knn_dists, rp_forest,) = nearest_neighbors(
        X,
        n_neighbors,
        "euclidean",
        {},
        False,
        random_state,
        verbose=False,
    )

    emb_graph, emb_sigmas, emb_rhos, emb_dists = fuzzy_simplicial_set(
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


@metric("density preservation", maximize=True)
def density_preservation(adata: AnnData) -> float:
    emb = adata.obsm["X_emb"]
    if np.any(np.isnan(emb)):
        return 0.0

    k, seed = 30, 42
    ro = _calculate_radii(
        adata.X.A if issparse(adata.X) else adata.X, n_neighbors=k, random_state=seed
    )
    re = _calculate_radii(emb, n_neighbors=k, random_state=seed)

    return pearsonr(ro, re)[0]
