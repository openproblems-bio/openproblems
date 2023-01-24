import anndata as ad
import scanpy as sc
from ivis import Ivis

# todo: allow using gpus instead!

## VIASH START
par = {
    "input": "resources_test/dimensionality_reduction/pancreas/train.h5ad",
    "output": "reduced.h5ad",
    "n_hvg": 1000,
    "n_pca_dims": 50
}
meta = {
    "functionality_name": "foo",
}
## VIASH END

print("Load input data", flush=True)
input = ad.read_h5ad(par["input"])
X_mat = input.layers["normalized"]

if par["n_hvg"]:
    print(f"Select top {par['n_hvg']} high variable genes", flush=True)
    idx = input.var["hvg_score"].to_numpy().argsort()[::-1][:par["n_hvg"]]
    X_mat = X_mat[:, idx]

print(f"Running PCA with {par['n_pca_dims']} dimensions", flush=True)
X_pca = sc.tl.pca(X_mat, n_comps=par["n_pca_dims"], svd_solver="arpack")[:, :2]

print("Run ivis")
# parameters taken from:
# https://bering-ivis.readthedocs.io/en/latest/scanpy_singlecell.html#reducing-dimensionality-using-ivis
ivis = Ivis(
    k=15,
    model="maaten",
    n_epochs_without_progress=5,
    verbose=0,
    embedding_dims=2,
)
X_emb = ivis.fit_transform(X_pca)

print("Create output AnnData", flush=True)
output = ad.AnnData(
    obs=input.obs[[]],
    obsm={
        "X_emb": X_emb
    },
    uns={
        "dataset_id": input.uns["dataset_id"],
        "normalization_id": input.uns["normalization_id"],
        "method_id": meta["functionality_name"]
    }
)

print("Write output to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")