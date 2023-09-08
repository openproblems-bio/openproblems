import anndata as ad
import numpy as np
import sklearn.decomposition
import sklearn.neighbors

## VIASH START
# The code between the the comments above and below gets stripped away before 
# execution. Here you can put anything that helps the prototyping of your script.
par = {
    "input_mod1": "resources_test/multimodal/integrated_mod1.h5ad",
    "input_mod2": "resources_test/multimodal/integrated_mod2.h5ad",
    "output": "resources_test/multimodal/score.h5ad",
    "proportion_neighbors": 0.1,
}

meta = {
    "functionality_name": "knn_auc"
}
## VIASH END

print("Reading adata file", flush=True)
adata_mod1 = ad.read_h5ad(par["input_mod1"])
adata_mod2 = ad.read_h5ad(par["input_mod2"])

print("Checking parameters", flush=True)
n_neighbors = int(np.ceil(par["proportion_neighbors"] * adata_mod1.layers["normalized"].shape[0]))

print("Compute KNN on PCA", flush=True)
_, indices_true = (
    sklearn.neighbors.NearestNeighbors(n_neighbors=n_neighbors)
    .fit(adata_mod1.obsm["X_svd"])
    .kneighbors(adata_mod1.obsm["X_svd"])
)

print("Compute KNN on integrated matrix", flush=True)
_, indices_pred = (
    sklearn.neighbors.NearestNeighbors(n_neighbors=n_neighbors)
    .fit(adata_mod1.obsm["integrated"])
    .kneighbors(adata_mod2.obsm["integrated"])
)

print("Check which neighbours match", flush=True)
print("Check which neighbours match", flush=True)
neighbors_match = np.zeros(n_neighbors, dtype=int)
for i in range(adata_mod1.layers["normalized"].shape[0]):
    _, pred_matches, true_matches = np.intersect1d(
        indices_pred[i], indices_true[i], return_indices=True
    )
    neighbors_match_idx = np.maximum(pred_matches, true_matches)
    neighbors_match += np.sum(
        np.arange(n_neighbors) >= neighbors_match_idx[:, None],
        axis=0,
    )

print("Compute area under neighbours match curve", flush=True)
print("Compute area under neighbours match curve", flush=True)
neighbors_match_curve = neighbors_match / (
    np.arange(1, n_neighbors + 1) * adata_mod1.layers["normalized"].shape[0]
)
area_under_curve = np.mean(neighbors_match_curve)

print("Store metic value", flush=True)
output_metric = ad.AnnData(
    layers={},
    obs=adata_mod1.obs[[]],
    var=adata_mod1.var[[]],
    uns={},
)

for key in adata_mod1.uns_keys():
    output_metric.uns[key] = adata_mod1.uns[key]

output_metric.uns["metric_ids"] = meta["functionality_name"]
output_metric.uns["metric_values"] = area_under_curve

print("Writing adata to file", flush=True)
output_metric.write_h5ad(par["output"], compression = "gzip")
