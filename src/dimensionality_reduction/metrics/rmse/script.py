import anndata as ad
from umap import UMAP, spectral
import scipy.spatial.distance as dist
from scipy.optimize import nnls
import numpy as np
from sklearn import decomposition, metrics

## VIASH START
par = {
    'input_reduced': 'reduced.h5ad',
    'input_test': 'test.h5ad',
    'output': 'score.h5ad',
}
meta = {
    'functionality_name': 'rmse',
}
## VIASH END

print("Load data", flush=True)
input_reduced = ad.read_h5ad(par['input_reduced'])
input_test = ad.read_h5ad(par['input_test'])

print('Reduce dimensionality of raw data', flush=True)
n_comps = 200
if not par['spectral']:
    input_reduced.obsm['high_dim'] = decomposition.TruncatedSVD(n_components = n_comps).fit_transform(input_test.layers['counts'])
    print('Compute RMSE between the full (or processed) data matrix and a dimensionally-reduced matrix, invariant to scalar multiplication', flush=True)
else:
    n_comps = min(n_comps, min(input_test.shape) - 2)
    graph = UMAP(transform_mode="graph").fit_transform(input_test.layers['counts'])
    input_reduced.obsm['high_dim'] = spectral.spectral_layout(
        input_test.layers['counts'], graph, n_comps, random_state=np.random.default_rng()
    )
    meta['functionality_name'] += ' spectral'
    print('Computes (RMSE) between high-dimensional Laplacian eigenmaps on the full (or processed) data matrix and the dimensionally-reduced matrix, invariant to scalar multiplication', flush=True)

high_dim_dist = dist.pdist(input_reduced.obsm['high_dim'])
low_dim_dist = dist.pdist(input_reduced.obsm["X_emb"])

scale, rmse = nnls(
        low_dim_dist[:, None], high_dim_dist
        )

print("Store metric value", flush=True)
input_reduced.uns['metric_ids'] =  meta['functionality_name']
input_reduced.uns['metric_values'] = rmse

print("Copy data to new AnnData object", flush=True)
output = ad.AnnData(
    uns={}
)
output.uns['normalization_id'] = input_reduced.uns['normalization_id']
output.uns['method_id'] = input_reduced.uns['method_id']
output.uns['dataset_id'] = input_reduced.uns['dataset_id']
output.uns['metric_ids'] =  input_reduced.uns['metric_ids']
output.uns['metric_values'] = input_reduced.uns['metric_values']

print("Write data to file", flush=True)
output.write_h5ad(par['output'], compression="gzip")