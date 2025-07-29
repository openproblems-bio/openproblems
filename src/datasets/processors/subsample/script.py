import scanpy as sc
import random
import numpy as np

### VIASH START
par = {
    "input": "resources_test/common/scicar_cell_lines/temp_mod1_full.h5ad",
    "input_mod2": "resources_test/common/scicar_cell_lines/temp_mod2_full.h5ad",
    "n_obs": 600,
    "n_vars": 1500,
    "keep_cell_type_categories": None,
    "keep_batch_categories": None,
    "keep_features": None,
    "keep_cell_type_categories": None,
    "keep_batch_categories": None,
    "even": False,
    "output": "subsample_mod1.h5ad",
    "output_mod2": "subsample_mod2.h5ad",
    "seed": 123
}
### VIASH END

if par["seed"]:
    print(f">> Setting seed to {par['seed']}", flush=True)
    random.seed(par["seed"])

print(">> Load data", flush=True)
adata_input = sc.read_h5ad(par["input"])

if par["input_mod2"] is not None:
    adata_mod2 = sc.read_h5ad(par["input_mod2"])

# copy counts to .X because otherwise filter_genes and filter_cells won't work
adata_input.X = adata_input.layers["counts"]
if par["input_mod2"] is not None:
    adata_mod2.X = adata_mod2.layers["counts"]

print(">> Determining output shape", flush=True)
min_obs_list = [par["n_obs"], adata_input.shape[0]]
if par["input_mod2"] is not None:
    min_obs_list.append(adata_mod2.shape[0])
n_obs = min(min_obs_list)

min_vars_list = [par["n_vars"], adata_input.shape[1]]
if par["input_mod2"] is not None:
    min_vars_list.append(adata_mod2.shape[1])
n_vars = min(min_vars_list)

print(">> Subsampling the observations", flush=True)
obs_filt = np.ones(dtype=np.bool_, shape=adata_input.n_obs)

# subset by cell_type
if par.get("keep_cell_type_categories"):
    print(f">> Selecting cell_type_categories {par['keep_cell_type_categories']}")
    obs_filt = obs_filt & adata_input.obs["cell_type"].isin(par["keep_cell_type_categories"])

# subset by batch
if par.get("keep_batch_categories"):
    print(f">> Selecting cell_type_categories {par['keep_batch_categories']}")
    obs_filt = obs_filt & adata_input.obs["batch"].isin(par["keep_batch_categories"])

# subsample evenly across batches or not
if par.get("even"):
    obs_evenly = "batch"
    choice_ix = np.where(obs_filt)[0]
    choice_batch = adata_input[choice_ix].obs[obs_evenly]
    names, counts = np.unique(choice_batch, return_counts=True)
    probs = dict(zip(names, 1 / counts / len(names)))
    
    choice_probs = [ probs[batch] for batch in choice_batch ]
    obs_index = np.random.choice(choice_ix, size=n_obs, replace=False, p=choice_probs)
else:
    obs_index = np.random.choice(np.where(obs_filt)[0], n_obs, replace=False)

# subsample obs
adata_output = adata_input[obs_index].copy()
if par["input_mod2"] is not None:
    adata_output_mod2 = adata_mod2[obs_index].copy()

# filter cells and genes
if par["input_mod2"] is not None:
    n_cells =  adata_output.X.sum(axis=1).A.flatten()
    n_cells_mod2 =  adata_output_mod2.X.sum(axis=1).A.flatten()
    keep_cells = np.minimum(n_cells, n_cells_mod2) > 1
    adata_output = adata_output[keep_cells, :].copy()
    adata_output_mod2 = adata_output_mod2[keep_cells, :].copy()

    sc.pp.filter_genes(adata_output, min_cells=1)
    sc.pp.filter_genes(adata_output_mod2, min_cells=1)
    
else:
    # todo: this should not remove features in keep_features!
    print(">> Remove empty observations and features", flush=True)
    sc.pp.filter_genes(adata_output, min_cells=1)
    sc.pp.filter_cells(adata_output, min_counts=2)

print(">> Subsampling the features", flush=True)
if par.get("keep_features"):
    initial_filt = adata_output.var_names.isin(par["keep_features"])
    initial_idx, *_ = initial_filt.nonzero()
    remaining_idx, *_ = (~initial_filt).nonzero()
    rest_idx = remaining_idx[np.random.choice(len(remaining_idx), n_vars - len(initial_idx), replace=False)]
    var_ix = np.concatenate([initial_idx, rest_idx])
else:
    var_ix = np.random.choice(adata_output.shape[1], n_vars, replace=False)
    if par["input_mod2"] is not None:
        var_ix_mod2 = np.random.choice(adata_output_mod2.shape[1], n_vars, replace=False)

#  subsample vars
adata_output = adata_output[:, var_ix].copy()
if par["input_mod2"] is not None:
    adata_output_mod2 = adata_output_mod2[:, var_ix_mod2].copy()

# filter cells and genes
if par["input_mod2"] is not None:
    n_cells =  adata_output.X.sum(axis=1).A.flatten()
    n_cells_mod2 =  adata_output_mod2.X.sum(axis=1).A.flatten()
    keep_cells = np.minimum(n_cells, n_cells_mod2) > 1
    adata_output = adata_output[keep_cells, :].copy()
    adata_output_mod2 = adata_output_mod2[keep_cells, :].copy()

    sc.pp.filter_genes(adata_output, min_cells=1)
    sc.pp.filter_genes(adata_output_mod2, min_cells=1)
    

else:
    # todo: this should not remove features in keep_features!
    print(">> Remove empty observations and features", flush=True)
    sc.pp.filter_genes(adata_output, min_cells=1)
    sc.pp.filter_cells(adata_output, min_counts=2)

print(">> Update dataset_id", flush=True)
adata_output.uns["dataset_id"] = adata_output.uns["dataset_id"] + "_subsample"
if par["input_mod2"] is not None:
    adata_output_mod2.uns["dataset_id"] = adata_output_mod2.uns["dataset_id"] + "_subsample"

# remove previously copied .X
del adata_output.X
if par["input_mod2"] is not None:
    del adata_output_mod2.X

print(">> Writing data", flush=True)
adata_output.write_h5ad(par["output"], compression=par["output_compression"])
if par["output_mod2"] is not None:
    adata_output_mod2.write_h5ad(par["output_mod2"], compression=par["output_compression"])
