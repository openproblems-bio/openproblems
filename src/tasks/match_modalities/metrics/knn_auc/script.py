import anndata as ad
import numpy as np
import sklearn.decomposition
import sklearn.neighbors

## VIASH START
par = {
  "input_integrated_mod1": "resources_test/match_modalities/scicar_cell_lines/integrated_mod1.h5ad",
  "input_integrated_mod2": "resources_test/match_modalities/scicar_cell_lines/integrated_mod2.h5ad",
  "input_solution_mod1": "resources_test/match_modalities/scicar_cell_lines/solution_mod1.h5ad",
  "input_solution_mod2": "resources_test/match_modalities/scicar_cell_lines/solution_mod2.h5ad",
  "output": "resources_test/multimodal/score.h5ad",
  "proportion_neighbors": 0.1,
}
meta = {
    "functionality_name": "knn_auc"
}
## VIASH END

print("Reading adata file", flush=True)
input_solution_mod1 = ad.read_h5ad(par["input_solution_mod1"])
input_solution_mod2 = ad.read_h5ad(par["input_solution_mod2"])

input_integrated_mod1 = ad.read_h5ad(par["input_integrated_mod1"])[input_solution_mod1.obs["permutation_indices"]]
input_integrated_mod2 = ad.read_h5ad(par["input_integrated_mod2"])[input_solution_mod2.obs["permutation_indices"]]

print("Checking parameters", flush=True)
n_neighbors = int(np.ceil(par["proportion_neighbors"] * input_solution_mod1.n_obs))

print("Compute KNN on PCA", flush=True)
_, indices_true = (
    sklearn.neighbors.NearestNeighbors(n_neighbors=n_neighbors)
    .fit(input_solution_mod1.obsm["X_svd"])
    .kneighbors(input_solution_mod2.obsm["X_svd"])
)

_, indices_pred = (
    sklearn.neighbors.NearestNeighbors(n_neighbors=n_neighbors)
    .fit(input_integrated_mod1.obsm["integrated"])
    .kneighbors(input_integrated_mod2.obsm["integrated"])
)

print("Check which neighbours match", flush=True)
neighbors_match = np.zeros(n_neighbors, dtype=int)
for i in range(input_solution_mod1.n_obs):
    _, pred_matches, true_matches = np.intersect1d(
        indices_pred[i], indices_true[i], return_indices=True
    )
    neighbors_match_idx = np.maximum(pred_matches, true_matches)
    neighbors_match += np.sum(
        np.arange(n_neighbors) >= neighbors_match_idx[:, None],
        axis=0,
    )

print("Compute area under neighbours match curve", flush=True)
neighbors_match_curve = neighbors_match / (
    np.arange(1, n_neighbors + 1) * input_solution_mod1.n_obs
)
area_under_curve = np.mean(neighbors_match_curve)

print("Store metric value", flush=True)
uns = {
  "dataset_id": input_solution_mod1.uns["dataset_id"],
  "normalization_id": input_solution_mod1.uns["normalization_id"],
  "method_id": input_integrated_mod1.uns["method_id"],
  "metric_ids": "knn_auc",
  "metric_values": area_under_curve
}
output_metric = ad.AnnData(
  shape=(0,0),
  uns=uns
)

print("Writing adata to file", flush=True)
output_metric.write_h5ad(par["output"], compression = "gzip")
