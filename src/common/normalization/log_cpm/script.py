import scanpy as sc

## VIASH START
par = {
    'input': "resources_test/label_projection/pancreas/dataset_subsampled.h5ad",
    'output': "output.h5ad"
}
meta = {
    "functionality_name": "normalize_log_cpm"
}
## VIASH END

print(">> Load data")
adata = sc.read_h5ad(par['input'])

print(">> Normalize data")
norm = sc.pp.normalize_total(adata, target_sum=1e6, key_added="size_factors", layer="counts", inplace=False)
lognorm = sc.pp.log1p(norm["X"])

print(">> Store output in adata")
adata.layers["lognorm"] = lognorm
adata.obs["norm_factor"] = norm["norm_factor"]
adata.uns["normalization_method"] = meta["functionality_name"].removeprefix("normalize_")

print(">> Write data")
adata.write_h5ad(par['output'], compression="gzip")
