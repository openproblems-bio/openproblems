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

with open(meta['config'], 'r', encoding="utf8") as file:
    config = yaml.safe_load(file)

output_type = config["functionality"]["info"]["output_type"]

print('Read input', flush=True)
input = ad.read_h5ad(par['input'])

if par['hvg']:
    print('Select HVGs', flush=True)
    input = input[:, input.var['hvg']].copy()

print('Run BBKNN', flush=True)
input.X = input.layers['normalized']
input = bbknn(input, batch='batch')
del input.X

print("Store outputs", flush=True)
input.uns['output_type'] = output_type
input.uns['hvg'] = par['hvg']
input.uns['method_id'] = meta['functionality_name']
input.write_h5ad(par['output'], compression='gzip')
