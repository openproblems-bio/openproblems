import scanpy as sc

## VIASH START
par = {
    'input': "resources_test/common/pancreas/dataset.h5ad",
    'output': "output.h5ad",
    'layer_output': "log_cp10k",
    'obs_size_factors': "log_cp10k_size_factors",
    'n_cp': 1e6,
}
meta = {
    "name": "normalize_log_cp10k"
}
## VIASH END

print(">> Load data", flush=True)
adata = sc.read_h5ad(par['input'])

print(">> Normalize data", flush=True)
if par["n_cp"] == -1:
    norm = sc.pp.normalize_total(
        adata, 
        target_sum=None, 
        layer="counts", 
        inplace=False
    )
else:
    norm = sc.pp.normalize_total(
        adata, 
        target_sum=par["n_cp"], 
        layer="counts", 
        inplace=False
    )
lognorm = sc.pp.log1p(norm["X"])

print(">> Store output in adata", flush=True)
adata.layers[par["layer_output"]] = lognorm
adata.obs[par["obs_size_factors"]] = norm["norm_factor"]
adata.uns["normalization_id"] = par["normalization_id"] or meta['name']

print(">> Write data", flush=True)
adata.write_h5ad(par['output'], compression="gzip")
