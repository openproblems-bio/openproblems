import anndata as ad
import numpy as np
from sklearn import manifold

## VIASH START
par = {
    "input_embedding": "resources_test/dimensionality_reduction/pancreas/reduced.h5ad",
    "input_solution": "resources_test/dimensionality_reduction/pancreas/test.h5ad",
    "output": "score.h5ad",
}
## VIASH END

print("Load data", flush=True)
input_solution = ad.read_h5ad(par["input_solution"])
input_embedding = ad.read_h5ad(par["input_embedding"])

high_dim = input_solution.layers["normalized"]
X_emb = input_embedding.obsm["X_emb"]

print("Reduce dimensionality of raw data", flush=True)
trustworthiness = manifold.trustworthiness(
    high_dim, X_emb, n_neighbors=15, metric="euclidean"
)

print("Create output AnnData object", flush=True)
output = ad.AnnData(
    uns={
        "dataset_id": input_solution.uns["dataset_id"],
        "normalization_id": input_solution.uns["normalization_id"],
        "method_id": input_embedding.uns["method_id"],
        "metric_ids": [ "trustworthiness" ],
        "metric_values": [ trustworthiness ]
    }
)

print("Write data to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")