import knn_smooth #NOTE: not a python package but just a py script on github
import numpy as np
import anndata as ad

## VIASH START
par = {
    'input_train': 'resources/label_projection/openproblems_v1/pancreas.split_dataset.output_train.h5ad',
    'output': 'output_knn.h5ad',
}
meta = {
    'functionality_name': 'foo',
}
## VIASH END

print("Load input data")
input_train = ad.read_h5ad(par["input_train"])

print("process data")
X = input_train.layers["counts"].transpose().toarray()
output_denoised = input_train.copy()
output_denoised.layers["counts"] = (knn_smooth.knn_smoothing(X, k=10)).transpose()

print("Writing data")
output_denoised.uns["method_id"] = par["functionality_name"]
output_denoised.write_h5ad(par["output"], compression="gzip")
