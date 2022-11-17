import anndata as ad

## VIASH START
par = {
    'input_train': 'output_train.h5ad',
    'input_test': 'output_test.h5ad',
    'output': 'output_baseline_PD.h5ad',
    'layer_input': 'counts',
    'baseline': 'perfect_denoising',
}
meta = {
    'functionality_name': 'foo',
}
## VIASH END

print("Load input data")
input_train = ad.read_h5ad(par['input_train'])
input_test = ad.read_h5ad(par['input_test'])


output_denoised = input_train.copy()

print("process data")
if par['baseline'] == "no_denoising":
    """Do nothing."""
    output_denoised.layers["denoised"] = input_train.layers[par['layer_input']].toarray()

elif par['baseline'] == "perfect_denoising":
    """Cheat."""
    output_denoised.layers["denoised"] = input_test.layers[par['layer_input']].toarray()

output_denoised.uns["method_id"] = meta['functionality_name']

print("Write Data")
output_denoised.write_h5ad(par['output'],compression="gzip")
