import anndata as ad
import numpy as np
from sklearn import manifold

## VIASH START
par = {
    'input_reduced': 'reduced.h5ad',
    'input_test': 'test.h5ad',
    'output': 'score.h5ad',
}
meta = {
    'functionality_name': 'trustworthiness',
}
## VIASH END

print("Load data", flush=True)
input_reduced = ad.read_h5ad(par['input_reduced'])
input_test = ad.read_h5ad(par['input_test'])

print('Reduce dimensionality of raw data', flush=True)
high_dim, low_dim = input_test.layers['counts'], input_reduced.obsm["X_emb"]
score = manifold.trustworthiness(
    high_dim, low_dim, n_neighbors=15, metric="euclidean"
)
# for large k close to #samples, it's higher than 1.0, e.g 1.0000073552559712
print("Store metric value", flush=True)
input_reduced.uns['metric_ids'] = meta['functionality_name']
input_reduced.uns['metric_values'] = float(np.clip(score, 0, 1))

print("Copy data to new AnnData object", flush=True)
output = ad.AnnData(
    uns={key: input_reduced.uns[key] for key in ["dataset_id", "normalization_id", "method_id", 'metric_ids', 'metric_values']}
)

print("Write data to file", flush=True)
output.write_h5ad(par['output'], compression="gzip")