import anndata as ad
import numpy as np
from scipy.sparse import issparse
from sklearn.decomposition import NMF

## VIASH START
par = {
  'input_single_cell': 'resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/single_cell_ref.h5ad',
  'input_spatial_masked': 'resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/spatial_masked.h5ad',
  'output': 'output.h5ad',
  'max_iter': 4000
}
meta = {
  'functionality_name': 'vanillanmf'
}
## VIASH END

print('Reading input files', flush=True)
input_single_cell = ad.read_h5ad(par['input_single_cell'])
input_spatial = ad.read_h5ad(par['input_spatial_masked'])

print('Generate predictions', flush=True)

n_types = input_single_cell.obs["cell_type"].cat.categories.shape[0]
vanila_nmf_model = NMF(
  n_components=n_types,
  beta_loss="kullback-leibler",
  solver="mu",
  max_iter=par['max_iter'],
  alpha_W=0.1,
  alpha_H=0.1,
  init="custom",
  random_state=42,
)

# Make profiles from single-cell expression dataset
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

X = input_spatial.layers['counts'].toarray()

Wa = vanila_nmf_model.fit_transform(
  X.astype(adata_means.X.dtype),
  H=adata_means.X,
  W=np.ones((input_spatial.shape[0], n_types), dtype=adata_means.X.dtype),
)

prop = Wa / Wa.sum(1)[:, np.newaxis]
input_spatial.obsm["proportions_pred"] = prop

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
