import anndata as ad
import scanpy as sc

## VIASH START
par = {
    "input": "resources_test/common/pancreas/train.h5ad",
    "output": "reduced.h5ad",
    "n_pca_dims": 50,
    "n_hvg": 1000
}
meta = {
    "functionality_name": "tsne",
    "config": "src/dimensionality_reduction/methods/tsne/config.vsh.yaml"
}
## VIASH END

print("Load input data", flush=True)
input = ad.read_h5ad(par["input"])

X_mat = input.layers["normalized"]

if par["n_hvg"]:
    print(f"Subsetting to {par['n_hvg']} HVG", flush=True)
    idx = input.var["hvg_score"].to_numpy().argsort()[::-1][:par["n_hvg"]]
    X_mat = X_mat[:, idx]

print("Computing PCA", flush=True)
input.obsm["X_pca"] = sc.tl.pca(X_mat, n_comps=par["n_pca_dims"], svd_solver="arpack")

print("Run t-SNE", flush=True)
sc.tl.tsne(input, use_rep="X_pca", n_pcs=par["n_pca_dims"])
X_emb = input.obsm["X_tsne"].copy()

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