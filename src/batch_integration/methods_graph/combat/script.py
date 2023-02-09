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

print("Run PCA", flush=True)
adata.obsm['X_emb'] = sc.pp.pca(
    adata.X,
    n_comps=50,
    use_highly_variable=False,
    svd_solver='arpack',
    return_info=False
)
del adata.X

print("Run KNN", flush=True)
sc.pp.neighbors(adata, use_rep='X_emb')

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['functionality_name']
adata.write_h5ad(par['output'], compression='gzip')
