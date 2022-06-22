##VIASH START
par = {
    'input': "../../resources/pancreas/toy_preprocessed_data.h5ad",
    'output': "output.h5ad"
}
##VIASH END

import scanpy as sc

print(">> Load data")
adata = sc.read(par['input'])

print(">> Normalize data")
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e6, key_added="size_factors")
sc.pp.log1p(adata)
adata.uns["normalization_method"] = "log_cpm"

print(">> Write data")
adata.write(par['output'], compression="gzip")
