import scanpy as sc
import random
import anndata as ad
import numpy as np

### VIASH START
par = {
    "input": "resources_test/common/pancreas/temp_dataset_full.h5ad",
    "keep_celltype_categories": None,
    "keep_batch_categories": None,
    "keep_features": ["HMGB2", "CDK1", "NUSAP1", "UBE2C"],
    "keep_celltype_categories": ["acinar", "beta"],
    "keep_batch_categories": ["celseq", "inDrop4", "smarter"],
    "even": True,
    "output": "toy_data.h5ad",
    "seed": 123
}
### VIASH END

if par["seed"]:
    print(f">> Setting seed to {par['seed']}", flush=True)
    random.seed(par["seed"])

print(">> Load data", flush=True)
adata_input = sc.read_h5ad(par["input"])

# copy counts to .X because otherwise filter_genes and filter_cells won't work
adata_input.X = adata_input.layers["counts"]

print(">> Determining output shape", flush=True)
n_obs = min(par["n_obs"], adata_input.shape[0])
n_vars = min(par["n_vars"], adata_input.shape[1])

print(">> Subsampling the observations", flush=True)
obs_filt = np.ones(dtype=np.bool_, shape=adata_input.n_obs)

# subset by celltype
if par.get("keep_celltype_categories"):
    print(f">> Selecting celltype_categories {par['keep_celltype_categories']}")
    obs_filt = obs_filt & adata_input.obs["celltype"].isin(par["keep_celltype_categories"])

# subset by batch
if par.get("keep_batch_categories"):
    print(f">> Selecting celltype_categories {par['keep_batch_categories']}")
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
    obs_index = np.random.choice(np.where(obs_filt)[0], n_vars, replace=False)

print(">> Subsampling the features", flush=True)
if par.get("keep_features"):
    initial_filt = adata_input.var_names.isin(par["keep_features"])
    initial_idx, *_ = initial_filt.nonzero()
    remaining_idx, *_ = (~initial_filt).nonzero()
    rest_idx = remaining_idx[np.random.choice(len(remaining_idx), n_vars - len(initial_idx), replace=False)]
    var_ix = np.concatenate([initial_idx, rest_idx])
else:
    var_ix = np.random.choice(adata_input.shape[1], n_vars, replace=False)


adata_output = adata_input[obs_index, var_ix].copy()

print(">> Remove empty observations and features", flush=True)
sc.pp.filter_genes(adata_output, min_cells=1)
sc.pp.filter_cells(adata_output, min_counts=2)

print(">> Update dataset_id", flush=True)
adata_output.uns["dataset_id"] = adata_output.uns["dataset_id"] + "_subsample"

# remove previously copied .X
del adata_output.X

print(">> Writing data")
adata_output.write_h5ad(par["output"])
