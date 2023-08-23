import yaml
import anndata as ad
from scib.integration import bbknn

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad',
    'hvg': True,
}
meta = {
    'functionality_name': 'foo',
    'config': 'bar'
}
## VIASH END

print('Read input', flush=True)
input = ad.read_h5ad(par['input'])

print('Run BBKNN', flush=True)
input.X = input.layers['normalized']
input = bbknn(input, batch='batch')
del input.X

print("Store outputs", flush=True)
input.uns['method_id'] = meta['functionality_name']
input.write_h5ad(par['output'], compression='gzip')
