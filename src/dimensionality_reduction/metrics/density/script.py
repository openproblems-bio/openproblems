import anndata as ad
import numpy as np
from typing import Optional
from umap.umap_ import fuzzy_simplicial_set
from umap.umap_ import nearest_neighbors
from umap import UMAP
from scipy.sparse import issparse
from scipy.stats import pearsonr

## VIASH START
par = {
    'input': 'output.h5ad',
    'output': 'score.h5ad',
}
meta = {
    'functionality_name': 'foo',
}
## VIASH END

print("Load data")
adata = ad.read_h5ad(par['input'])

print('Reduce dimensionality of raw data')
_K = 30  # number of neighbors
_SEED = 42
_, ro, _ = UMAP(
    n_neighbors=_K, random_state=_SEED, densmap=True, output_dens=True
).fit_transform(adata.layers['counts'])
# in principle, we could just call _calculate_radii(high_dim, ...)
# this is made sure that the test pass (otherwise, there was .02 difference in corr)
(knn_indices, knn_dists, rp_forest,) = nearest_neighbors(
    adata.obsm['X_emb'],
    _K,
    "euclidean",
    {},
    False,
    _SEED,
    verbose=False,
)

emb_graph, emb_sigmas, emb_rhos, emb_dists = fuzzy_simplicial_set(
    adata.obsm['X_emb'],
    _K,
    _SEED,
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
re = np.log(epsilon + (re / mu_sum))

adata.uns['metric_ids'] = 'density'
adata.uns['metric_values'] = pearsonr(ro, re)[0]

print("Write data to file")
adata.write_h5ad(par['output'], compression="gzip")