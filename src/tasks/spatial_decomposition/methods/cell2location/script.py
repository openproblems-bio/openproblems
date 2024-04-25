import anndata as ad
import numpy as np
from cell2location.cluster_averages.cluster_averages import compute_cluster_averages
from cell2location.models import Cell2location
from cell2location.models import RegressionModel
from torch.nn import ELU

## VIASH START
par = {
  'input_single_cell': 'resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/single_cell_ref.h5ad',
  'input_spatial_masked': 'resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/spatial_masked.h5ad',
  'output': 'output.h5ad',
  'detection_alpha': 20.0,
  'n_cells_per_location': 20,
  'hard_coded_reference': True,
  'amortised': False,
  'num_samples': 1000,
  'sc_batch_size': 2500,
  'st_batch_size': None,
  'max_epochs_sc': 250,
  'max_epochs_st': 5000
}
meta = {
  'functionality_name': 'cell2location'
}
## VIASH END

print('Reading input files', flush=True)
input_single_cell = ad.read_h5ad(par['input_single_cell'])
input_spatial = ad.read_h5ad(par['input_spatial_masked'])

input_single_cell.X = input_single_cell.layers["counts"]
input_spatial.X = input_spatial.layers["counts"]

if not par["hard_coded_reference"]:
  if "batch" in input_single_cell.obs.columns:
      input_single_cell.obs["batch_key"] = input_single_cell.obs["batch"].copy()
  else:
    input_single_cell.obs["batch_key"] = "all"
  # REFERENCE SIGNATURE ESTIMATION FROM scRNA
  # prepare anndata for the regression model
  RegressionModel.setup_anndata(
    adata=input_single_cell,
    # 10X reaction / sample / batch
    batch_key="batch_key",
    # cell type, covariate used for constructing signatures
    labels_key="cell_type",
  )
  sc_model = RegressionModel(input_single_cell)
  sc_model.train(max_epochs=par["max_epochs_sc"], batch_size=par["sc_batch_size"])
  # In this section, we export the estimated cell abundance
  # (summary of the posterior distribution).
  input_single_cell = sc_model.export_posterior(
    input_single_cell,
    sample_kwargs={"num_samples": par["num_samples"], "batch_size": par["sc_batch_size"]},
  )
  # export estimated expression in each cluster
  try:
    means_per_cluster = input_single_cell.varm["means_per_cluster_mu_fg"]
  except KeyError:
    # sometimes varm fails for unknown reason
    means_per_cluster = input_single_cell.var
  means_per_cluster = means_per_cluster[
    [
      f"means_per_cluster_mu_fg_{i}"
      for i in input_single_cell.uns["mod"]["factor_names"]
    ]
  ].copy()
  means_per_cluster.columns = input_single_cell.uns["mod"]["factor_names"]
else:
  means_per_cluster = compute_cluster_averages(
    input_single_cell,
    labels="cell_type",
    layer=None,
    use_raw=False,
  )

# SPATIAL MAPPING
# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(input_spatial.var_names, means_per_cluster.index)
input_spatial = input_spatial[:, intersect].copy()
means_per_cluster = means_per_cluster.loc[intersect, :].copy()

# prepare anndata for cell2location model
input_spatial.obs["sample"] = "all"
Cell2location.setup_anndata(adata=input_spatial, batch_key="sample")
cell2location_kwargs = dict(
    cell_state_df=means_per_cluster,
    # the expected average cell abundance: tissue-dependent hyper-prior which can be estimated from paired histology:
    # here = average in the simulated dataset
    N_cells_per_location=par["n_cells_per_location"],
    # hyperparameter controlling normalisation of within-experiment variation in RNA detection:
    detection_alpha=par["detection_alpha"],
)
if par["amortised"]:
    cell2location_kwargs["amortised"] = True
    cell2location_kwargs["encoder_mode"] = "multiple"
    cell2location_kwargs["encoder_kwargs"] = {
        "dropout_rate": 0.1,
        "n_hidden": {
            "single": 256,
            "n_s_cells_per_location": 10,
            "b_s_groups_per_location": 10,
            "z_sr_groups_factors": 64,
            "w_sf": 256,
            "detection_y_s": 20,
        },
        "use_batch_norm": False,
        "use_layer_norm": True,
        "n_layers": 1,
        "activation_fn": ELU,
    }
# create and train the model
st_model = Cell2location(input_spatial, **cell2location_kwargs)
st_model.train(
    max_epochs=par["max_epochs_st"],
    # train using full data (batch_size=None)
    batch_size=par["st_batch_size"],
    # use all data points in training because we need to estimate cell abundance at all locations
    train_size=1,
)
# In this section, we export the estimated cell abundance (summary of the posterior distribution).
input_spatial = st_model.export_posterior(
    input_spatial,
    sample_kwargs={
        "num_samples": par["num_samples"],
        "batch_size": par["st_batch_size"],
    },
)

input_spatial.obsm["proportions_pred"] = input_spatial.obsm["q05_cell_abundance_w_sf"].values
input_spatial.obsm["proportions_pred"] /= input_spatial.obsm["proportions_pred"].sum(axis=1)[:, None]

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

print("Write output to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")