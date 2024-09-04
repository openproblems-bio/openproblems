import anndata as ad
from scvi.external import RNAStereoscope
from scvi.external import SpatialStereoscope

## VIASH START
par = {
  'input_single_cell': 'resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/single_cell_ref.h5ad',
  'input_spatial_masked': 'resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/spatial_masked.h5ad',
  'output': 'output.h5ad', 
  'max_epochs_sc': 100,
  'max_epochs_sp': 1000
}
meta = {
  'functionality_name': 'stereoscope'
}
## VIASH END

print('Reading input files', flush=True)
input_single_cell = ad.read_h5ad(par['input_single_cell'])
input_spatial = ad.read_h5ad(par['input_spatial_masked'])

input_single_cell.X = input_single_cell.layers["counts"]
input_spatial.X = input_spatial.layers["counts"]

print('Generate predictions', flush=True)

RNAStereoscope.setup_anndata(input_single_cell, labels_key="cell_type")
sc_model = RNAStereoscope(input_single_cell)
sc_model.train(
  max_epochs=par["max_epochs_sc"],
  # early_stopping=True,
  # early_stopping_monitor="elbo_validation"
)

SpatialStereoscope.setup_anndata(input_spatial)
st_model = SpatialStereoscope.from_rna_model(input_spatial, sc_model)
st_model.train(
  max_epochs=par["max_epochs_sp"],
  # early_stopping=True,
  # early_stopping_monitor="elbo_validation"
)
input_spatial.obsm["proportions_pred"] = st_model.get_proportions().to_numpy()

print("Write output AnnData to file", flush=True)
output = ad.AnnData(
  obs=input_spatial.obs[[]],
  var=input_spatial.var[[]],
  obsm={
    'coordinates': input_spatial.obsm['coordinates'],
    'proportions_pred': input_spatial.obsm['proportions_pred']
  },
  layers={
    'counts': input_spatial.layers['counts']
  },
  uns={
    'cell_type_names': input_spatial.uns['cell_type_names'],
    'dataset_id': input_spatial.uns['dataset_id'],
    'method_id': meta['functionality_name']
  }
)
output.write_h5ad(par['output'], compression='gzip')
