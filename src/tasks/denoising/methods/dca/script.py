import anndata as ad
from dca.api import dca
import scipy

## VIASH START
par = {
    'input_train': 'resources_test/denoising/pancreas/train.h5ad',
    'output': 'output_dca.h5ad',
    'epochs': 300,
}
meta = {
    'functionality_name': 'dca',
}
## VIASH END

print("load input data", flush=True)
input_train = ad.read_h5ad(par['input_train'])

print("move layer to X", flush=True)
input_dca = ad.AnnData(X=input_train.layers["counts"])
del input_train.X

print("running dca", flush=True)
dca(input_dca, epochs=par["epochs"])

print("moving X back to layer", flush=True)
input_train.layers["denoised"] = scipy.sparse.csr_matrix(input_dca.X)

print("Writing data", flush=True)
input_train.uns["method_id"] = meta["functionality_name"]
input_train.write_h5ad(par["output"], compression="gzip")
