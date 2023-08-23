import yaml
import anndata as ad
from scib.integration import scanorama

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
adata = ad.read_h5ad(par['input'])

print('Run scanorama', flush=True)
adata.X = adata.layers['normalized']
adata.obsm['X_emb'] = scanorama(adata, batch='batch').obsm['X_emb']
del adata.X

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['functionality_name']
adata.write(par['output'], compression='gzip')
