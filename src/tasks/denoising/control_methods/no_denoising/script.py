import anndata as ad

## VIASH START
par = {
    'input_train': 'output_train.h5ad',
    'output': 'output_ND.h5ad',
}
meta = {
    'functionality_name': 'foo',
}
## VIASH END

print("Load input data", flush=True)
input_train = ad.read_h5ad(par['input_train'])

print("Process data", flush=True)
input_train.layers["denoised"] = input_train.layers['counts']

input_train.uns["method_id"] = meta['functionality_name']

print("Write Data", flush=True)
input_train.write_h5ad(par['output'],compression="gzip")
