import anndata as ad
from dca.api import dca

# NOTE: pickup later. DCA package uses old keras imports ?

## VIASH START
par = {
    'input_train': 'output_train.h5ad',
    'output': 'output_dca.h5ad',
    'epochs': 300,
}
meta = {
    'functionality_name': 'dca',
}
## VIASH END

print("load input data")
input_train = ad.read_h5ad(par['input_train'])

print("process data")

# make adata object with train counts
# run DCA
dca(input_train, epochs=par["epochs"])

# set denoised to Xmat
output_denoised = input_train.copy
output_denoised.X = input_train.X

print("Writing data")
output_denoised.uns["method_id"] = meta['functionality_name']
output_denoised.write_h5ad(par['output'], compression="gzip")




