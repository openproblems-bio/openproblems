import anndata as ad
import numpy as np
from scipy import sparse

## VIASH START
par = {
  "input_integrated_mod1": "resources_test/match_modalities/scicar_cell_lines/integrated_mod1.h5ad",
  "input_integrated_mod2": "resources_test/match_modalities/scicar_cell_lines/integrated_mod2.h5ad",
  "input_solution_mod1": "resources_test/match_modalities/scicar_cell_lines/solution_mod1.h5ad",
  "input_solution_mod2": "resources_test/match_modalities/scicar_cell_lines/solution_mod2.h5ad",
  "output": "resources_test/multimodal/score.h5ad",
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

print("Computing MSE", flush=True)
def _square(X):
	if sparse.issparse(X):
		X.data = X.data ** 2
		return X
	else:
		return X ** 2


X = input_integrated_mod1.obsm["integrated"].toarray()
Y = input_integrated_mod2.obsm["integrated"].toarray()

X_shuffled = X[np.random.permutation(np.arange(X.shape[0])), :]
error_random = np.mean(np.sum(_square(X_shuffled - Y)))
error_abs = np.mean(np.sum(_square(X - Y)))
metric_value = (error_abs / error_random).item()

print("Store metric value", flush=True)
uns = {
  "dataset_id": input_solution_mod1.uns["dataset_id"],
  "normalization_id": input_solution_mod1.uns["normalization_id"],
  "method_id": input_integrated_mod1.uns["method_id"],
  "metric_ids": "mse",
  "metric_values": metric_value
}
output_metric = ad.AnnData(
  shape=(0,0),
  uns=uns
)

print("Writing adata to file", flush=True)
output_metric.write_h5ad(par["output"], compression = "gzip")
