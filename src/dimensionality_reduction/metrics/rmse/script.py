import anndata as ad
import scipy.spatial.distance as dist
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

print("Load data")
input_reduced = ad.read_h5ad(par['input_reduced'])
input_test = ad.read_h5ad(par['input_test'])

print('Reduce dimensionality of raw data')
input_reduced.obsm['svd'] = decomposition.TruncatedSVD(n_components = 200).fit_transform(input_test.layers['counts'])

print('Compute pairwise distance between points in a matrix and format it into a squared-form vector.')
high_dim_dist_matrix = dist.squareform(dist.pdist(input_reduced.obsm['svd']))
low_dim_dist_matrix = dist.squareform(dist.pdist(input_reduced.obsm["X_emb"]))

print('Compute RMSE between the full (or processed) data matrix and a dimensionally-reduced matrix')
y_actual = high_dim_dist_matrix
y_predict = low_dim_dist_matrix
rmse = np.sqrt(metrics.mean_squared_error(y_actual, y_predict))

print('Compute Kruskal stress between the full (or processed) data matrix and a dimensionally-reduced matrix')
diff = high_dim_dist_matrix - low_dim_dist_matrix
kruskal_matrix = np.sqrt(diff**2 / sum(low_dim_dist_matrix**2))
kruskal_score = np.sqrt(sum(diff**2) / sum(low_dim_dist_matrix**2))

print("Store metric value")
input_reduced.uns['metric_ids'] =  meta['functionality_name']
input_reduced.uns['metric_values'] = rmse

print("Delete obs matrix")
del input_reduced.obsm

print("Write data to file")
input_reduced.write_h5ad(par['output'], compression="gzip")