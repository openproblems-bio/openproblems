# TODO: this should be a output_type: features method.

import scanpy as sc
from scipy.sparse import csr_matrix

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad',
    'hvg': True
}

meta = {
    'funcionality_name': 'foo'
}

## VIASH END

print('Read input', flush=True)
adata = sc.read_h5ad(par['input'])

if par['hvg']:
    print('Select HVGs', flush=True)
    adata = adata[:, adata.var['hvg']].copy()

print('Run Combat', flush=True)
adata.X = adata.layers['normalized']
adata.X = sc.pp.combat(adata, key='batch', inplace=False)
adata.X = csr_matrix(adata.X)

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['functionality_name']
adata.write_h5ad(par['output'], compression='gzip')
