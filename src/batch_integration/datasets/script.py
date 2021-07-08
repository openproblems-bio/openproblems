## VIASH START
par = {
    'adata': '../resources/pancreas.h5ad',
    'name': 'pancreas',
    'label': 'celltype',
    'batch': 'tech',
    'hvgs': 2000,
    'output': 'adata_out.h5ad',
    'debug': True
}
resources_dir = '../../common/utils/'
## VIASH END

print('Importing libraries')
import scanpy as sc
import sys
import os

sys.path.append(os.path.abspath(resources_dir))

from preprocessing import log_scran_pooling

if par['debug']:
    import pprint

    pprint.pprint(par)

adata_file = par['adata']
name = par['name']
label = par['label']
batch = par['batch']
hvgs = par['hvgs']
output = par['output']

print('Read adata')
adata = sc.read(adata_file)
assert name == adata.uns['name']

# Rename columns
adata.obs['labels'] = adata.obs[label]
adata.obs['batch'] = adata.obs[batch]
adata.layers['counts'] = adata.X

print(f'Select {hvgs} highly variable genes')
if adata.n_obs > hvgs:
    sc.pp.subsample(adata, n_obs=hvgs)

print('Normalisation with scran')
log_scran_pooling(adata)
adata.layers['logcounts'] = adata.X

print('Transformation: PCA')
sc.tl.pca(
    adata,
    svd_solver='arpack',
    return_info=True,
)
adata.obsm['X_uni'] = adata.obsm['X_pca']

print('Transformation: kNN')
sc.pp.neighbors(adata, use_rep='X_uni', key_added='uni')

print('Writing adata to file')
adata.write(output, compression='gzip')
