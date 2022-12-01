import anndata as ad
import numpy as np
from sklearn import manifold

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
high_dim, low_dim = adata.layers['counts'], adata.obsm["X_emb"]
score = manifold.trustworthiness(
    high_dim, low_dim, n_neighbors=15, metric="euclidean"
)
# for large k close to #samples, it's higher than 1.0, e.g 1.0000073552559712
adata.uns['metric_ids'] = 'trustworthiness'
adata.uns['metric_values'] = float(np.clip(score, 0, 1))

print("Write data to file")
adata.write_h5ad(par['output'], compression="gzip")