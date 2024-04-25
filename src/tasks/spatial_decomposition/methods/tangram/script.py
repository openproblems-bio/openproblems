import anndata as ad
import pandas as pd
import scanpy as sc
import tangram as tg
import torch

## VIASH START
par = {
  'input_single_cell': 'resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/single_cell_ref.h5ad',
  'input_spatial_masked': 'resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/spatial_masked.h5ad',
  'output': 'output.h5ad',
  'num_epochs': 1000,
  'n_markers': 100
}
meta = {
  'functionality_name': 'tangram'
}
## VIASH END

print('Reading input files', flush=True)
input_single_cell = ad.read_h5ad(par['input_single_cell'])
input_spatial = ad.read_h5ad(par['input_spatial_masked'])

print('Generate predictions', flush=True)
# analysis based on github.com/broadinstitute/Tangram/blob/master/tutorial_tangram_with_squidpy.ipynb
# using tangram from PyPi, not github version

input_single_cell.X = input_single_cell.layers["counts"]
input_spatial.X = input_spatial.layers["counts"]

# pre-process single cell data
sc.pp.normalize_total(input_single_cell, 1e4)
sc.pp.log1p(input_single_cell)
# identify marker genes
sc.tl.rank_genes_groups(input_single_cell, groupby="cell_type", use_raw=False)

# extract marker genes to data frame
markers_df = pd.DataFrame(input_single_cell.uns["rank_genes_groups"]["names"]).iloc[0:par['n_markers'], :]

# get union of all marker genes
markers = list(set(markers_df.melt().value.values))

# match genes between single cell and spatial data
tg.pp_adatas(input_single_cell, input_spatial, genes=markers)

# get device
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# map single cells to spatial locations
ad_map = tg.map_cells_to_space(
  input_single_cell,
  input_spatial,
  device=device,
  num_epochs=par['num_epochs'],
)

# transfer labels from mapped cells to spatial location
tg.project_cell_annotations(adata_map=ad_map, adata_sp=input_spatial, annotation="cell_type")

# normalize scores
pred_props = input_spatial.obsm["tangram_ct_pred"].to_numpy()
input_spatial.obsm["proportions_pred"] = pred_props / pred_props.sum(axis=1)[:, None]

# remove un-normalized predictions
del input_spatial.obsm["tangram_ct_pred"]

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
