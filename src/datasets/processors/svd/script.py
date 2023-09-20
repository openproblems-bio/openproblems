import anndata as ad
import sklearn.decomposition


## VIASH START
par = {
    "input": "resources_test/common/scicar_cell_lines/normalized_mod1.h5ad",
    "input_mod2": "resources_test/common/scicar_cell_lines/normalized_mod2.h5ad",
    "output": "output.h5ad",
    "layer_input": "normalized",
    "obsm_embedding": "X_svd",
    "num_components": 100,
}
## VIASH END

print(">> Load data", flush=True)
adata = ad.read(par["input"])
if par["input_mod2"] is not None:
    adata2 = ad.read(par["input_mod2"])

print(">> check parameters", flush=True)
min_list = [par["num_components"], min(adata.layers[par["layer_input"]].shape) - 1]

if par["input_mod2"] is not None:
    min_list.append(min(adata2.layers[par["layer_input"]].shape) - 1)

n_svd = min(min_list)


print(">> Run SVD", flush=True)
svd1 = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata.layers[par["layer_input"]])
if par["input_mod2"] is not None:
    svd2 = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata2.layers[par["layer_input"]])

print(">> Storing output", flush=True)
adata.obsm[par["obsm_embedding"]] = svd1
if par["input_mod2"] is not None:
    adata2.obsm[par["obsm_embedding"]] = svd2


print(">> Writing data", flush=True)
adata.write_h5ad(par["output"])
if par["input_mod2"] is not None:
    adata2.write_h5ad(par["output_mod2"])

