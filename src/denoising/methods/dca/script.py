import anndata as ad
from dca.api import dca

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

print("move layer to X")
input_train.X = input_train.layers["counts"]
del input_train.layers["counts"]

print("running dca")
dca(input_train, epochs=par["epochs"])

print("moving X back to layer")
input_train.layers["denoised"] = input_train.X
del input_train.X

print("Writing data")
input_train.uns["method_id"] = meta["functionality_name"]
input_train.write_h5ad(par["output"], compression="gzip")
