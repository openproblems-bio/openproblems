import yaml
import anndata as ad
from scib.integration import scanvi

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad',
    'hvg': True,
}
meta = {
    'functionality_name' : 'foo',
    'config': 'bar'
}
## VIASH END



print('Read input', flush=True)
adata = ad.read_h5ad(par['input'])

if par['hvg']:
    print('Select HVGs', flush=True)
    adata = adata[:, adata.var['hvg']].copy()

print('Run scanvi', flush=True)
adata.X = adata.layers['normalized']
adata = scanvi(adata, batch='batch', labels='label')
del adata.X

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['functionality_name']
adata.write_h5ad(par['output'], compression='gzip')
