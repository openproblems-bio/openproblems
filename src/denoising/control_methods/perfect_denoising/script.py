import anndata as ad

## VIASH START
par = {
    'input_train': 'resources_test/denoising/pancreas/train.h5ad',
    'input_test': 'resources_test/denoising/pancreas/test.h5ad',
    'output': 'output_baseline_PD.h5ad',
}
meta = {
    'functionality_name': 'foo',
}
## VIASH END

print("Load input data")
input_train = ad.read_h5ad(par['input_train'])
input_test = ad.read_h5ad(par['input_test'])

print("Process data")
input_train.layers["denoised"] = input_test.layers['counts']

input_train.uns["method_id"] = meta['functionality_name']

print("Write Data")
input_train.write_h5ad(par['output'],compression="gzip")
