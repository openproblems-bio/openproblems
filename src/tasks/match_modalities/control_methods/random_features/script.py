import anndata as ad
import numpy as np

## VIASH START

par = {
    "input_mod1": "resources_test/common/scicar_cell_lines/dataset_mod1.h5ad",
    "input_mod2": "resources_test/common/scicar_cell_lines/dataset_mod2.h5ad",
    "output_mod1": "output.mod1.h5ad",
    "output_mod2": "output.mod2.h5ad",
}

meta = {
    "functionality_name": "random_features"
}

## VIASH END

print("Reading input h5ad file", flush=True)
adata_mod1 = ad.read_h5ad(par["input_mod1"])
adata_mod2 = ad.read_h5ad(par["input_mod2"])

print("Generating random features", flush=True)
adata_mod1.obsm["integrated"] = adata_mod1.obsm["X_svd"][np.random.permutation(np.arange(adata_mod1.shape[0]))]
adata_mod2.obsm["integrated"] = adata_mod1.obsm["X_svd"][np.random.permutation(np.arange(adata_mod1.shape[0]))]

print("Write output to file", flush=True)
adata_mod1.uns["method_id"] = meta["functionality_name"]
adata_mod2.uns["method_id"] = meta["functionality_name"]
adata_mod1.write_h5ad(par["output_mod1"], compression="gzip")
adata_mod2.write_h5ad(par["output_mod2"], compression="gzip")