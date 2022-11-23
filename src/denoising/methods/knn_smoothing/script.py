import knn_smooth
import numpy as np
import anndata as ad

## VIASH START
par = {
    'input_train': 'output/denoising/pancreas.split_data.output_train.h5ad',
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
input_train.layers["denoised"] = (knn_smooth.knn_smoothing(X, k=10)).transpose()

print("Writing data")
input_train.uns["method_id"] = meta["functionality_name"]
input_train.write_h5ad(par["output"], compression="gzip")
