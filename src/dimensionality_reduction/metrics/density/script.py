import anndata as ad
import numpy as np
from typing import Optional
from umap import UMAP
from scipy.sparse import issparse
from scipy.stats import pearsonr

## VIASH START
par = {
    'input_reduced': 'reduced.h5ad',
    'input_test': 'test.h5ad',
    'output': 'score.h5ad',
}
meta = {
    'functionality_name': 'density',
}
## VIASH END
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

print("Load data")
input_reduced = ad.read_h5ad(par['input_reduced'])
input_test = ad.read_h5ad(par['input_test'])

if np.any(np.isnan(input_reduced.obsm['X_emb'])):
    density = 0.0
else:
    _K = 30  # number of neighbors
    _SEED = 42
    
    print('Reduce dimensionality of raw data', flush=True)
    _, ro, _ = UMAP(
        n_neighbors=_K, random_state=_SEED, densmap=True, output_dens=True
    ).fit_transform(input_test.layers['counts'])
    re = _calculate_radii(input_reduced.obsm['X_emb'], n_neighbors=_K, random_state=_SEED)
    
    print('Compute Density between the full (or processed) data matrix and a dimensionally-reduced matrix', flush=True)
    density = pearsonr(ro, re)[0]

print("Store metric value", flush=True)
input_reduced.uns['metric_ids'] =  meta['functionality_name']
input_reduced.uns['metric_values'] = density

print("Copy data to new AnnData object", flush=True)
output = ad.AnnData(
    uns={key: input_reduced.uns[key] for key in ["dataset_id", "normalization_id", "method_id", 'metric_ids', 'metric_values']}
)

print("Write data to file", flush=True)
output.write_h5ad(par['output'], compression="gzip")