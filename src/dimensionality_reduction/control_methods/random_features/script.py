import anndata as ad
import numpy as np

## VIASH START
par = {
    "input": "resources_test/dimensionality_reduction/pancreas/test.h5ad",
    "output": "reduced.h5ad",
}
meta = {
    "functionality_name": "random_features",
}
## VIASH END

print("Load input data", flush=True)
input = ad.read_h5ad(par["input"])

print("Create random embedding", flush=True)
X_emb = np.random.normal(0, 1, (input.shape[0], 2))

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