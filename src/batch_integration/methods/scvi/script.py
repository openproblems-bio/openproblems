import yaml
import anndata as ad
from scib.integration import scvi

## VIASH START
par = {
    'input': 'resources_test/batch_integration/unintegrated.h5ad',
    'output': 'output.h5ad',
    'hvg': True,
}
meta = {
    'functionality_name' : 'foo',
    'config': 'bar'
}
## VIASH END

with open(meta['config'], 'r', encoding="utf8") as file:
    config = yaml.safe_load(file)

output_type = config["functionality"]["info"]["output_type"]

print('Read input', flush=True)
adata = ad.read_h5ad(par['input'])

if par['hvg']:
    print('Select HVGs', flush=True)
    adata = adata[:, adata.var['hvg']].copy()

print('Run scvi', flush=True)
adata.X = adata.layers['normalized']
adata = scvi(adata, batch='batch')
del adata.X

print("Store outputs", flush=True)
adata.uns['output_type'] = output_type
adata.uns['hvg'] = par['hvg']
adata.uns['method_id'] = meta['functionality_name']
adata.write_h5ad(par['output'], compression='gzip')
