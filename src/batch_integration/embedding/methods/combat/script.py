import scanpy as sc
from scipy.sparse import csr_matrix

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/processed.h5ad',
    'output': 'output.h5ad',
    'hvg': True,
    'scaling': True
}

meta = {
    'functionality_name' : 'foo'
}
## VIASH END

adata_file = par['input']
output = par['output']
hvg = par['hvg']
scaling = par['scaling']

print('Read input', flush=True)
adata = sc.read_h5ad(adata_file)

if hvg:
    print('Select HVGs', flush=True)
    adata = adata[:, adata.var['highly_variable']].copy()

if scaling:
    print('Scale', flush=True)
    adata.X = adata.layers['logcounts_scaled']
else:
    adata.X = adata.layers['logcounts']

print('Integrate', flush=True)
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
    obs= adata.obs[[]],
    obsm={
        'X_emb': X_emb
    },
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': meta['functionality_name'],
        'scaled': par['scaling'],
        'hvg': par('hvg')
    },
    layers = adata.layers
)

print('Write to output', flush=True)
adata.write_h5ad(par['output'], compression='gzip')
