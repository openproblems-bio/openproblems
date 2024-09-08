import anndata as ad

## VIASH START
par = {
  'input_single_cell': 'resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/single_cell.h5ad',
  'input_spatial_masked': 'resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/spatial_masked.h5ad',
  'input_solution': 'resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/solution.h5ad',
  'output': 'output.h5ad'
}
meta = {
  'functionality_name': 'true_proportions'
}
## VIASH END

print('Reading input files', flush=True)
input_single_cell = ad.read_h5ad(par['input_single_cell'])
input_spatial_masked = ad.read_h5ad(par['input_spatial_masked'])
input_solution = ad.read_h5ad(par['input_solution'])

print('Generate predictions', flush=True)
input_spatial_masked.obsm["proportions_pred"] = input_solution.obsm["proportions_true"]

print("Write output AnnData to file", flush=True)
output = ad.AnnData(
  obs=input_spatial_masked.obs[[]],
  var=input_spatial_masked.var[[]],
  uns={
    'cell_type_names': input_spatial_masked.uns['cell_type_names'],
    'dataset_id': input_spatial_masked.uns['dataset_id'],
    'method_id': meta['functionality_name']
  },
  obsm={
    'coordinates': input_spatial_masked.obsm['coordinates'],
    'proportions_pred': input_spatial_masked.obsm['proportions_pred']
  },
  layers={
    'counts': input_spatial_masked.layers['counts']
  }
)
output.write_h5ad(par['output'], compression='gzip')
