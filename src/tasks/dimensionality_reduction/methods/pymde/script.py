import anndata as ad
import scanpy as sc
import pymde

## VIASH START
par = {
    "input": "resources_test/dimensionality_reduction/pancreas/dataset.h5ad",
    "output": "reduced.h5ad",
    "embed_method": "neighbors",
    "n_hvg": 1000,
    "n_pca_dims": 50,
}
meta = {
    "functionality_name": "foo",
}
## VIASH END

if par["embed_method"] == "neighbors":
    mde_fn = pymde.preserve_neighbors
elif par["embed_method"] == "distances":
    mde_fn = pymde.preserve_distances
else:
    raise ValueError(f"Unknown embedding method: {par['embed_method']}")

print("Load input data", flush=True)
input = ad.read_h5ad(par["input"])
X_mat = input.layers["normalized"]

if par["n_hvg"]:
    print(f"Select top {par['n_hvg']} high variable genes", flush=True)
    idx = input.var["hvg_score"].to_numpy().argsort()[::-1][:par["n_hvg"]]
    X_mat = X_mat[:, idx]

print(f"Compute PCA", flush=True)
X_pca = sc.tl.pca(X_mat, n_comps=par["n_pca_dims"], svd_solver="arpack")

print(f"Run MDE", flush=True)
X_emb = (
    mde_fn(X_pca, embedding_dim=2, verbose=True)
    .embed(verbose=True)
    .detach()
    .numpy()
)

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