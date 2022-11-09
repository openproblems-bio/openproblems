## VIASH START
# The code between the the comments above and below gets stripped away before 
# execution. Here you can put anything that helps the prototyping of your script.
par = {
    "input": "out_bash/modality_alignment/methods/citeseq_cbmc_mnn.h5ad",
    "output": "out_bash/modality_alignment/metrics/citeseq_cbmc_mnn_knn_auc.h5ad"
}
## VIASH END

print("Importing libraries")
import anndata
import scprep
import numpy as np
from scipy import sparse

print("Reading adata file")
adata = anndata.read_h5ad(par["input"])

print("Computing MSE")
def _square(X):
	if sparse.issparse(X):
		X.data = X.data ** 2
		return X
	else:
		return scprep.utils.toarray(X) ** 2

X = scprep.utils.toarray(adata.obsm["aligned"])
Y = scprep.utils.toarray(adata.obsm["mode2_aligned"])

X_shuffled = X[np.random.permutation(np.arange(X.shape[0])), :]
error_random = np.mean(np.sum(_square(X_shuffled - Y)))
error_abs = np.mean(np.sum(_square(X - Y)))
metric_value = error_abs / error_random

print("Store metic value")
adata.uns["metric_id"] = "mse"
adata.uns["metric_value"] = metric_value

print("Writing adata to file")
adata.write_h5ad(par["output"], compression = "gzip")
