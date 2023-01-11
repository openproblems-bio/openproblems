## VIASH START
par = {
    'input': './src/batch_integration/resources/datasets_pancreas.h5ad',
    'output': './src/batch_integration/resources/pancreas_bbknn.h5ad',
    'hvg': True,
    'scaling': True,
    'debug': True
}
## VIASH END

print('Importing libraries')
from pprint import pprint
import scanpy as sc
from scipy.sparse import csr_matrix

if par['debug']:
    pprint(par)

adata_file = par['input']
output = par['output']
hvg = par['hvg']
scaling = par['scaling']

print('Read adata')
adata = sc.read_h5ad(adata_file)

if hvg:
    print('Select HVGs')
    adata = adata[:, adata.var['highly_variable']].copy()

if scaling:
    print('Scale')
    adata.X = adata.layers['logcounts_scaled']
else:
    adata.X = adata.layers['logcounts']

print('Integrate')
adata.X = sc.pp.combat(adata, key='batch', inplace=False)
adata.X = csr_matrix(adata.X)

print('Postprocess data')
adata.obsm['X_emb'] = sc.pp.pca(
    adata.X,
    n_comps=50,
    use_highly_variable=False,
    svd_solver='arpack',
    return_info=False
)
sc.pp.neighbors(adata, use_rep='X_emb')

print('Save HDF5')
adata.uns['hvg'] = hvg
adata.uns['scaled'] = scaling

adata.write(output, compression='gzip')
