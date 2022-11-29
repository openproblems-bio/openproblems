import anndata as ad
import scipy.spatial.distance as dist
import numpy as np
from sklearn import decomposition, metrics

## VIASH START
par = {
    'input': 'resources_test/common/pancreas/dataset.h5ad',
    'output': 'output.h5ad',
}
meta = {
    'functionality_name': 'foo',
}
## VIASH END

print("Load data")
adata = ad.read_h5ad(par['input'])

print('Reduce dimensionality of raw data')
adata.obsm['svd'] = decomposition.TruncatedSVD(n_components = 200).fit_transform(adata.layers['counts'])

print('Compute pairwise distance between points in a matrix and format it into a squared-form vector.')
high_dim_dist_matrix = dist.squareform(dist.pdist(adata.obsm['svd']))
low_dim_dist_matrix = dist.squareform(dist.pdist(adata.obsm["X_emb"]))

print('Compute RMSE between the full (or processed) data matrix and a dimensionally-reduced matrix')
y_actual = high_dim_dist_matrix
y_predict = low_dim_dist_matrix
rmse = np.sqrt(metrics.mean_squared_error(y_actual, y_predict))

print('Compute Kruskal stress between the full (or processed) data matrix and a dimensionally-reduced matrix')
diff = high_dim_dist_matrix - low_dim_dist_matrix
kruskal_matrix = np.sqrt(diff**2 / sum(low_dim_dist_matrix**2))
kruskal_score = np.sqrt(sum(diff**2) / sum(low_dim_dist_matrix**2))

print("Store metric value")
if 'metric_ids' not in adata.uns.keys():
    adata.uns['metric_ids'] = []
    adata.uns['metric_values'] = []

adata.uns['metric_ids'] += ['kruskal', 'rmse']
adata.uns['metric_values'] += [kruskal_score, rmse]
adata.obsm['kruskal'] = kruskal_matrix

print("Write data to file")
adata.write_h5ad(par['output'], compression="gzip")