import anndata as ad
from dca.api import dca
from scipy import sparse

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

print("load input data")
input_train = ad.read_h5ad(par['input_train'])

print("running dca")
mod = dca(input_train, epochs=par["epochs"], copy=True)

print("Writing data")
input_train.layers["denoised"] = mod.layers["counts"]

print("Writing data")
input_train.uns["method_id"] = meta["functionality_name"]
input_train.write_h5ad(par["output"], compression="gzip")
