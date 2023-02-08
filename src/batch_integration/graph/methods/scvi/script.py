# TODO: this should be a output_type: embedding method.
import scanpy as sc
from scib.integration import scvi

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad',
    'hvg': True,
}
meta = {
    'functionality_name' : 'foo'
}
## VIASH END

print('Read input', flush=True)
adata = sc.read_h5ad(par['input'])

if par['hvg']:
    print('Select HVGs', flush=True)
    adata = adata[:, adata.var['hvg']].copy()

print('Run scvi', flush=True)
adata.X = adata.layers['normalized']
adata = scvi(adata, batch='batch')

print('Run kNN', flush=True)
sc.pp.neighbors(adata, use_rep='X_emb')

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['functionality_name']
adata.write_h5ad(par['output'], compression='gzip')
