import anndata as ad
from muon import prot as pt

## VIASH START
par = {
    'input': "resources_test/common/openproblems_neurips2021/bmmc_cite/dataset_mod2.h5ad",
    'output': "output_norm.h5ad"
}
meta = {
    'functionality_name': "clr"
}
## VIASH END

print("Load data", flush=True)
adata = ad.read_h5ad(par['input'])

print("Normalize data", flush=True)
input_adata = ad.AnnData(X=adata.layers["counts"])
normalized_counts = pt.pp.clr(input_adata, inplace=False)
if not normalized_counts:
    raise RuntimeError("CLR failed to return the requested output layer")

print("Store output in adata", flush=True)
adata.layers[par["layer_output"]] = normalized_counts.X
adata.uns["normalization_id"] = par["normalization_id"] or meta['functionality_name']

print("Write data", flush=True)
adata.write_h5ad(par['output'], compression="gzip")
