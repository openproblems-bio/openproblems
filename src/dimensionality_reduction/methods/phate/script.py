import anndata as ad
from phate import PHATE

## VIASH START
par = {
    "input": "resources_test/dimensionality_reduction/pancreas/train.h5ad",
    "output": "reduced.h5ad",
    "n_pca_dims": 50,
    "n_hvg": 1000,
    "gamma": 1
}
meta = {
    "functionality_name": "foo",
}
## VIASH END

print("Load input data", flush=True)
input = ad.read_h5ad(par["input"])

X_mat = input.layers["normalized"]

if par["n_hvg"]:
    print(f"Subsetting to {par['n_hvg']} HVG", flush=True)
    idx = input.var["hvg_score"].to_numpy().argsort()[::-1][:par["n_hvg"]]
    X_mat = X_mat[:, idx]

print("Run PHATE", flush=True)
phate_op = PHATE(n_pca=par["n_pca_dims"], verbose=False, n_jobs=-1, gamma=par["gamma"])
X_emb = phate_op.fit_transform(X_mat)

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