import anndata as ad
import numpy as np
from scipy.optimize import nnls
from scipy.sparse import issparse

## VIASH START
par = {
  'input_single_cell': 'resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/single_cell_ref.h5ad',
  'input_spatial_masked': 'resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/spatial_masked.h5ad',
  'output': 'output.h5ad'
}
meta = {
  'functionality_name': 'nnls'
}
## VIASH END

print('Reading input files', flush=True)
input_single_cell = ad.read_h5ad(par['input_single_cell'])
input_spatial = ad.read_h5ad(par['input_spatial_masked'])

# Compute means over each 'cell_type'
labels = input_single_cell.obs['cell_type'].cat.categories
n_var = input_single_cell.shape[1]
means = np.empty((labels.shape[0], n_var))
for i, lab in enumerate(labels):
  adata_lab = input_single_cell[input_single_cell.obs['cell_type'] == lab]
  x_lab = adata_lab.layers['counts']
  means[i, :] = x_lab.mean(axis=0).flatten()
adata_means = ad.AnnData(means)
adata_means.obs_names = labels
adata_means.var_names = input_single_cell.var_names

X = adata_means.X.T
y = input_spatial.layers['counts'].T
res = np.zeros((y.shape[1], X.shape[1]))  # (voxels, cells)
for i in range(y.shape[1]):
  x, _ = nnls(X, y[:, i].toarray().reshape(-1))
  res[i] = x

# Normalize coefficients to sum to 1
res[res < 0] = 0
res = res / res.sum(axis=1, keepdims=1)

input_spatial.obsm["proportions_pred"] = res

print("Write output AnnData to file", flush=True)
output = ad.AnnData(
  obs=input_spatial.obs[[]],
  var=input_spatial.var[[]],
  uns={
    'cell_type_names': input_spatial.uns['cell_type_names'],
    'dataset_id': input_spatial.uns['dataset_id'],
    'method_id': meta['functionality_name']
  },
  obsm={
    'coordinates': input_spatial.obsm['coordinates'],
    'proportions_pred': input_spatial.obsm['proportions_pred']
  },
  layers={
    'counts': input_spatial.layers['counts']
  }
)
output.write_h5ad(par['output'], compression='gzip')
