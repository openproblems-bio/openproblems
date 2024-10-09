import anndata as ad
from muon import atac as ac

## VIASH START
par = {
    'input': "resources_test/common/openproblems_neurips2021/bmmc_cite/dataset_mod2.h5ad",
    'output': "output_norm.h5ad"
}
meta = {
    'name': "tfidf"
}
## VIASH END

print("Load data", flush=True)
adata = ad.read_h5ad(par['input'])

print("Normalize data", flush=True)
input_adata = ad.AnnData(X=adata.layers["counts"])
normalized_counts = ac.pp.tfidf(input_adata, inplace=False)

print("Store output in adata", flush=True)
adata.layers[par["layer_output"]] = normalized_counts
adata.uns["normalization_id"] = par["normalization_id"] or meta['name']

print("Write data", flush=True)
adata.write_h5ad(par['output'], compression="gzip")
