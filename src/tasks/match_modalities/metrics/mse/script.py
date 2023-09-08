import anndata as ad
import scprep
import numpy as np
from scipy import sparse

## VIASH START
# The code between the the comments above and below gets stripped away before 
# execution. Here you can put anything that helps the prototyping of your script.
par = {
    "input_mod1": "resources_test/multimodal/integrated_mod1.h5ad",
    "input_mod2": "resources_test/multimodal/integrated_mod2.h5ad",
    "output": "resources_test/multimodal/mse.h5ad"
}
## VIASH END

print("Reading adata file", flush=True)
adata_mod1 = ad.read_h5ad(par["input_mod1"])
adata_mod2 = ad.read_h5ad(par["input_mod2"])

print("Computing MSE", flush=True)
def _square(X):
	if sparse.issparse(X):
		X.data = X.data ** 2
		return X
	else:
		return scprep.utils.toarray(X) ** 2

X = scprep.utils.toarray(adata_mod1.obsm["integrated"])
Y = scprep.utils.toarray(adata_mod2.obsm["integrated"])

X_shuffled = X[np.random.permutation(np.arange(X.shape[0])), :]
error_random = np.mean(np.sum(_square(X_shuffled - Y)))
error_abs = np.mean(np.sum(_square(X - Y)))
metric_value = error_abs / error_random

output_metric = ad.AnnData(
	layers={},
	obs=adata_mod1.obs[[]],
	var=adata_mod1.var[[]],
	uns={}
)

for key in adata_mod1.uns_keys():
    output_metric.uns[key] = adata_mod1.uns[key]

print("Store metic value", flush=True)
output_metric.uns["metric_ids"] = meta["functionality_name"]
output_metric.uns["metric_values"] = metric_value

print("Writing adata to file", flush=True)
output_metric.write_h5ad(par["output"], compression = "gzip")
