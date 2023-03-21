from typing import Any, Callable, Dict, Tuple
import openproblems as op
import scanpy as sc
import scipy
import pandas as pd

## VIASH START
par = {
    "dataset_id": "scicar_mouse_kidney",
    "obs_celltype": "celltype",
    "obs_batch": "replicate",
    "obs_tissue": None,
    "layer_counts": "counts",
    "output": "test_data.h5ad",
}
meta = {
    "resources_dir": "src/datasets/loaders/openproblems_v1/"
}
## VIASH END


# make dataset lookup table
# If need be, this could be stored in a separate yaml file
dataset_funs: Dict[str, Tuple[Callable, Dict[str, Any]]] = {
    "citeseq_cbmc": (op.data.multimodal.citeseq.load_citeseq_cbmc, {}),
    "scicar_cell_lines": (op.data.multimodal.scicar.load_scicar_cell_lines, {}),
    "scicar_mouse_kidney": (op.data.multimodal.scicar.load_scicar_mouse_kidney, {}),
}

# fetch dataset
dataset_fun, kwargs = dataset_funs[par["dataset_id"]]

print("Fetch dataset", flush=True)
adata = dataset_fun(**kwargs)

print(f"source adata: {adata}", flush=True)

# construct modality2 dataset
mod2_var_data = {
    key.replace("mode2_var_", ""): adata.uns[key]
    for key in adata.uns.keys()
    if key.startswith("mode2_var_")
}
mod2_var = pd.DataFrame(
    mod2_var_data,
    index=adata.uns["mode2_var"]
)
mod2_obs = adata.obs.loc[adata.uns["mode2_obs"]]
mod2 = sc.AnnData(
    obs=mod2_obs,
    var=mod2_var,
    layers={ "counts": adata.obsm["mode2"] }
)

# construct modality1 dataset
mod1 = adata.copy()
mod1.uns = { key: value for key, value in mod1.uns.items() if not key.startswith("mode2_")}
mod1.obsm = { key: value for key, value in mod1.obsm.items() if not key.startswith("mode2_")}
mod1.obsp = { key: value for key, value in mod1.obsp.items() if not key.startswith("mode2_")}
mod1.varm = { key: value for key, value in mod1.varm.items() if not key.startswith("mode2_")}
mod1.varp = { key: value for key, value in mod1.varp.items() if not key.startswith("mode2_")}

# override values one by one because adata.uns and
# metadata are two different classes.
for key, value in dataset_fun.metadata.items():
    print(f"Setting .uns['{key}']", flush=True)
    mod1.uns[key] = value
    mod2.uns[key] = value

print("Setting .obs['celltype']", flush=True)
if par["obs_celltype"]:
    if par["obs_celltype"] in mod1.obs:
        mod1.obs["celltype"] = mod1.obs[par["obs_celltype"]]
        mod2.obs["celltype"] = mod2.obs[par["obs_celltype"]]
    else:
        print(f"Warning: key '{par['obs_celltype']}' could not be found in adata.obs.", flush=True)

print("Setting .obs['batch']", flush=True)
if par["obs_batch"]:
    if par["obs_batch"] in mod1.obs:
        mod1.obs["batch"] = mod1.obs[par["obs_batch"]]
        mod2.obs["batch"] = mod2.obs[par["obs_batch"]]
    else:
        print(f"Warning: key '{par['obs_batch']}' could not be found in adata.obs.", flush=True)

print("Setting .obs['tissue']", flush=True)
if par["obs_tissue"]:
    if par["obs_tissue"] in mod1.obs:
        mod1.obs["tissue"] = mod1.obs[par["obs_tissue"]]
        mod2.obs["tissue"] = mod2.obs[par["obs_tissue"]]
    else:
        print(f"Warning: key '{par['obs_tissue']}' could not be found in adata.obs.", flush=True)

if par["layer_counts"] and par["layer_counts"] in mod1.layers:
    print(f"Temporarily moving mod1.layers['{par['layer_counts']}']", flush=True)
    mod1_X = mod1.layers[par["layer_counts"]]
    del mod1.layers[par["layer_counts"]]
else:
    print("Temporarily moving mod1.X", flush=True)
    mod1_X = mod1.X
    del mod1.X

if par["sparse"] and not scipy.sparse.issparse(mod1_X):
    print("Make mod1 counts sparse", flush=True)
    mod1_X = scipy.sparse.csr_matrix(mod1_X)

if par["sparse"] and not scipy.sparse.issparse(mod2.layers["counts"]):
    print("Make mod2 counts sparse", flush=True)
    mod2.layers["counts"] = scipy.sparse.csr_matrix(mod2.layers["counts"])

print("Moving .X to .layers['counts']", flush=True)
mod1.layers["counts"] = mod1_X

# just in case
del mod1.X
del mod2.X

print("Add metadata to uns", flush=True)
metadata_fields = [
    "dataset_id", "dataset_name", "data_url", "data_reference",
    "dataset_summary", "dataset_description" "dataset_organism"
]
uns_metadata = {
    id: par[id]
    for id in metadata_fields
    if id in par
}
mod1.uns.update(uns_metadata)
mod2.uns.update(uns_metadata)

print("Writing adata to file", flush=True)
mod1.write_h5ad(par["output_mod1"], compression="gzip")
mod2.write_h5ad(par["output_mod2"], compression="gzip")
