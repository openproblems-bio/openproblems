# TODO: this should be a output_type: feature method.
import scanpy as sc
from scib.integration import scanorama

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad',
    'hvg': True,
}
meta = {
    'functionality_name': 'foo'
}
## VIASH END

print('Read input', flush=True)
adata = sc.read_h5ad(par['input'])

if par['hvg']:
    print('Select HVGs', flush=True)
    adata = adata[:, adata.var['hvg']]

print('Run scanorama', flush=True)
adata.X = adata.layers['normalized']
adata.X = scanorama(adata, batch='batch').X

print("Run PCA", flush=True)
sc.pp.pca(
    adata,
    n_comps=50,
    use_highly_variable=False,
    svd_solver='arpack',
    return_info=True
)
del adata.X

print("Run KNN", flush=True)
sc.pp.neighbors(adata, use_rep='X_pca')

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['functionality_name']
adata.write_h5ad(par['output'], compression='gzip')
