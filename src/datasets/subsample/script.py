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
    # "keep_celltype_categories": ["acinar", "beta"],
    # "keep_batch_categories": ["celseq", "inDrop4", "smarter"],
    "even": True,
    "output": "toy_data.h5ad",
    "seed": 123
}
### VIASH END

def filter_genes_cells(adata):
    """Remove empty cells and genes."""
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.filter_cells(adata, min_counts=2)


def subsample_even(adata, n_obs, even_obs):
    """Subsample a dataset evenly across an obs.
    Parameters
    ----------
    adata : AnnData
    n_obs : int
        Total number of cells to retain
    even_obs : str
        `adata.obs[even_obs]` to be subsampled evenly across partitions.
    Returns
    -------
    adata : AnnData
        Subsampled AnnData object
    """
    import scanpy as sc

    values = adata.obs[even_obs].unique()
    adatas = []
    n_obs_per_value = n_obs // len(values)
    for v in values:
        adata_subset = adata[adata.obs[even_obs] == v].copy()
        sc.pp.subsample(adata_subset, n_obs=min(n_obs_per_value, adata_subset.shape[0]))
        adatas.append(adata_subset)

    adata_out = ad.concat(adatas, label="_obs_batch")

    adata_out.uns = adata.uns
    adata_out.varm = adata.varm
    adata_out.varp = adata.varp
    return adata_out


if par["seed"]:
    print(f">> Setting seed to {par['seed']}", flush=True)
    random.seed(par["seed"])

print(">> Load data", flush=True)
adata_input = sc.read_h5ad(par["input"])

# copy counts to .X because otherwise filter_genes and filter_cells won't work
adata_input.X = adata_input.layers["counts"]

# filter by celltype
if par.get("keep_celltype_categories"):
    print(f">> Selecting celltype_categories {par['keep_celltype_categories']}")
    idx = adata_input.obs["celltype"].isin(par["keep_celltype_categories"])
    adata_input = adata_input[idx]

# filter by batch
if par.get("keep_batch_categories"):
    print(f">> Selecting celltype_categories {par['keep_batch_categories']}")
    idx = adata_input.obs["batch"].isin(par["keep_batch_categories"])
    adata_input = adata_input[idx]

print(">> Remove empty observations and features", flush=True)
filter_genes_cells(adata_input)

print(">> Subsampling the observations", flush=True)
n_obs = min(500, adata_input.shape[0])
n_vars = min(500, adata_input.shape[1])
if par.get("even"):
    adata_output = subsample_even(adata_input, n_obs, "batch")
else:
    adata_output = sc.pp.subsample(adata_input, n_obs=n_obs, copy=True)


print(">> Subsampling the features", flush=True)
if par.get("keep_features"):
    initial_filt = adata_output.var_names.isin(par["keep_features"])
    initial_idx, *_ = initial_filt.nonzero()
    remaining_idx, *_ = (~initial_filt).nonzero()
    rest_idx = remaining_idx[np.random.choice(len(remaining_idx), n_vars - len(initial_idx), replace=False)]
    feature_ix = np.concatenate([initial_idx, rest_idx])
    adata_output = adata_output[:, feature_ix]
else:
    feature_ix = np.random.choice(adata_input.shape[1], n_vars, replace=False)
    adata_output = adata_output[:, feature_ix]

print(">> Remove empty observations and features", flush=True)
filter_genes_cells(adata_output)

print(">> Update dataset_id", flush=True)
adata_output.uns["dataset_id"] = adata_output.uns["dataset_id"] + "_subsample"

# remove previously copied .X
del adata_output.X

print(">> Writing data")
adata_output.write_h5ad(par["output"])
