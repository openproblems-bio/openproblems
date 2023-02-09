import scanpy as sc
from scipy.sparse import csr_matrix

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/processed.h5ad',
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
    adata = adata[:, adata.var['highly_variable']].copy()

print('Run Combat', flush=True)
adata.X = adata.layers['normalized']
adata.X = sc.pp.combat(adata, key='batch', inplace=False)
adata.X = csr_matrix(adata.X)

print('Postprocess data', flush=True)
X_emb = sc.pp.pca(
    adata.X,
    n_comps=50,
    use_highly_variable=False,
    svd_solver='arpack',
    return_info=False
)

print('Create output AnnData object', flush=True)
output = sc.AnnData(
    obs= adata.obs,
    obsm={
        'X_emb': X_emb
    },
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': meta['functionality_name'],
        'hvg': par['hvg']
    },
    layers = {key: value for key, value in adata.layers.items()}
)

print('Write to output', flush=True)
output.write_h5ad(par['output'], compression='gzip')
