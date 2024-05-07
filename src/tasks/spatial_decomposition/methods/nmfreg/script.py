import anndata as ad
import numpy as np
from scipy.optimize import nnls
from scipy.sparse import issparse
from sklearn.decomposition import NMF
from sklearn.preprocessing import StandardScaler

## VIASH START
par = {
  'input_single_cell': 'resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/single_cell_ref.h5ad',
  'input_spatial_masked': 'resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/spatial_masked.h5ad',
  'output': 'output.h5ad',
  'n_components': 30
}
meta = {
  'functionality_name': 'nmfreg'
}
## VIASH END

print('Reading input files', flush=True)
input_single_cell = ad.read_h5ad(par['input_single_cell'])
input_spatial = ad.read_h5ad(par['input_spatial_masked'])

n_types = input_single_cell.obs["cell_type"].cat.categories.shape[0]

# Learn from reference
X = input_single_cell.layers['counts']
X_norm = X / X.sum(1)
X_scaled = StandardScaler(with_mean=False).fit_transform(X_norm)
model = NMF(
  n_components=par['n_components'],
  init="random",
  random_state=42
)
Ha = model.fit_transform(X_scaled)
Wa = model.components_

cluster_df = input_single_cell.obs[["cell_type"]].copy()
cluster_df.loc[:, "factor"] = np.argmax(Ha, axis=1)
cluster_df.loc[:, "code"] = cluster_df.cell_type.values.codes
factor_to_cluster_map = np.array(
  [
    np.histogram(
      cluster_df.loc[cluster_df.factor == k, "code"],
      bins=n_types,
      range=(0, n_types),
    )[0]
    for k in range(par['n_components'])
  ]
).T

factor_to_best_celltype = np.argmax(factor_to_cluster_map, axis=0)

factor_to_best_celltype_matrix = np.zeros((par['n_components'], n_types))
for i, j in enumerate(factor_to_best_celltype):
  factor_to_best_celltype_matrix[i, j] = 1

Ha_norm = StandardScaler(with_mean=False).fit_transform(Ha)
sc_deconv = np.dot(Ha_norm**2, factor_to_best_celltype_matrix)
sc_deconv = sc_deconv / sc_deconv.sum(1)[:, np.newaxis]

# Start run on actual spatial data
X_sp = input_spatial.layers['counts']
X_sp_norm = X_sp / X_sp.sum(1)
X_sp_scaled = StandardScaler(with_mean=False).fit_transform(X_sp_norm)

bead_prop_soln = np.array([nnls(Wa.T, X_sp_scaled[b, : ].toarray().reshape(-1))[0] for b in range(X_sp_scaled.shape[0])])
bead_prop_soln = StandardScaler(with_mean=False).fit_transform(bead_prop_soln)
bead_prop = np.dot(bead_prop_soln, factor_to_best_celltype_matrix)

prop = bead_prop / bead_prop.sum(1)[:, np.newaxis]
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