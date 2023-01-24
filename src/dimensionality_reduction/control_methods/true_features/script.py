import anndata as ad

## VIASH START
par = {
    "input": "resources_test/dimensionality_reduction/pancreas/test.h5ad",
    "output": "reduced.h5ad",
    "n_hvg": 100,
    "use_normalized_layer": False
}
meta = {
    "functionality_name": "true_features",
}
## VIASH END

print("Load input data", flush=True)
input = ad.read_h5ad(par["input"])

print("Create high dimensionally embedding with all features", flush=True)
if par["use_normalized_layer"]:
    X_emb = input.layers["counts"].toarray()
else:
    X_emb = input.layers["normalized"].toarray()

if par["n_hvg"]:
    print(f"Select top {par['n_hvg']} high variable genes", flush=True)
    idx = input.var["hvg_score"].to_numpy().argsort()[::-1][:par["n_hvg"]]
    X_emb = X_emb[:, idx]

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