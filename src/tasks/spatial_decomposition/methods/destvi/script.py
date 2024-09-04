import anndata as ad
from scvi.model import CondSCVI
from scvi.model import DestVI

## VIASH START
par = {
  'input_single_cell': 'resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/single_cell_ref.h5ad',
  'input_spatial_masked': 'resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/spatial_masked.h5ad',
  'output': 'output.h5ad', 
  'max_epochs_sc': 500,
  'max_epochs_sp': 5000
}
meta = {
  'functionality_name': 'destvi'
}
## VIASH END

print('Reading input files', flush=True)
input_single_cell = ad.read_h5ad(par['input_single_cell'])
input_spatial = ad.read_h5ad(par['input_spatial_masked'])

input_single_cell.X = input_single_cell.layers["counts"]
input_spatial.X = input_spatial.layers["counts"]

CondSCVI.setup_anndata(input_single_cell, labels_key="cell_type")
sc_model = CondSCVI(input_single_cell, weight_obs=False)
sc_model.train(
  max_epochs=par['max_epochs_sc'],
  early_stopping=True,
  train_size=0.9,
  validation_size=0.1,
  early_stopping_monitor="elbo_validation",
)

DestVI.setup_anndata(input_spatial)
st_model = DestVI.from_rna_model(input_spatial, sc_model)
st_model.train(
  max_epochs=par['max_epochs_sp'],
  batch_size=min(int(input_spatial.n_obs / 20 + 3), 128),
  plan_kwargs={"min_kl_weight": 3.0, "max_kl_weight": 3},
)
input_spatial.obsm["proportions_pred"] = st_model.get_proportions().to_numpy()

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

print("Write output to file", flush=True)
output.write_h5ad(par['output'], compression='gzip')