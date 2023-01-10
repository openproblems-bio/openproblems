import anndata as ad
import scanpy as sc

## VIASH START
par = {
    "input": "resources_test/dimensionality_reduction/pancreas/train.h5ad",
    "output": "reduced.h5ad",
    "n_hvg": 1000,
}
meta = {
    "functionality_name": "pca",
}
## VIASH END

print("Load input data", flush=True)
input = ad.read_h5ad(par["input"])
X_mat = input.layers["normalized"]

if par["n_hvg"]:
    print(f"Select top {par['n_hvg']} high variable genes", flush=True)
    idx = input.var["hvg_score"].to_numpy().argsort()[::-1][:par["n_hvg"]]
    X_mat = X_mat[:, idx]

print(f"Running PCA", flush=True)
X_emb = sc.tl.pca(X_mat, n_comps=2, svd_solver="arpack")[:, :2]

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