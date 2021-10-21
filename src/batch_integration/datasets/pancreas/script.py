## VIASH START
import os

print(os.getcwd())
par = {
    'adata': './src/batch_integration/datasets/resources/data_loader_pancreas.h5ad',
    'label': 'celltype',
    'batch': 'tech',
    'hvgs': 2000,
    'output': './src/batch_integration/datasets/resources/datasets_pancreas.h5ad',
    'debug': True
}
resources_dir = './src/batch_integration/datasets'
## VIASH END

print('Importing libraries')
import scanpy as sc
# import sys
# sys.path.append(resources_dir)
# from utils import log_scran_pooling

if par['debug']:
    import pprint

    pprint.pprint(par)

adata_file = par['adata']
label = par['label']
batch = par['batch']
hvgs = par['hvgs']
output = par['output']

print('Read adata')
adata = sc.read(adata_file)

# Rename columns
adata.obs.rename(columns={label: 'label', batch: 'batch'}, inplace=True)
adata.layers['counts'] = adata.X

print(f'Select {hvgs} highly variable genes')
if adata.n_obs > hvgs:
    sc.pp.subsample(adata, n_obs=hvgs)

#print('Normalisation with scran')
#log_scran_pooling(adata)
#adata.layers['logcounts'] = adata.X

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
