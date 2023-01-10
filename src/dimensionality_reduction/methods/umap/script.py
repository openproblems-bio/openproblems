import anndata as ad
from umap import UMAP
import scanpy as sc

## VIASH START
par = {
    "input": "resources_test/common/pancreas/train.h5ad",
    "output": "reduced.h5ad",
    "n_pca_dims": False,
}
meta = {
    "functionality_name": "densmap",
    "config": "src/dimensionality_reduction/methods/densmap/config.vsh.yaml"
}
## VIASH END

print("Load input data", flush=True)
input = ad.read_h5ad(par["input"])
X_mat = input.layers["normalized"]

if par["n_hvg"]:
    print(f"Select top {par['n_hvg']} high variable genes", flush=True)
    idx = input.var["hvg_score"].to_numpy().argsort()[::-1][:par["n_hvg"]]
    X_mat = X_mat[:, idx]

if par["n_pca_dims"]:
    print("Apply PCA to normalized data", flush=True)
    umap_input = sc.tl.pca(
        X_mat,
        n_comps=par["n_pca_dims"],
        svd_solver="arpack"
    )
else:
    print("Use normalized data as input for UMAP", flush=True)
    umap_input = X_mat

print("Run UMAP", flush=True)
X_emb = UMAP(densmap=False, random_state=42).fit_transform(umap_input)

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