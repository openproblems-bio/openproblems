## VIASH START
par = {
    'adata': './src/batch_integration/resources/datasets_pancreas.h5ad',
    'output': './src/batch_integration/resources/pancreas_bbknn.h5ad',
    'hvg': 100,
    'scaling': True,
    'debug': True
}
## VIASH END

print('Importing libraries')
import pprint
import scanpy as sc
from scIB.integration import runBBKNN
from scIB.preprocessing import hvg_batch
from scIB.preprocessing import scale_batch

if par['debug']:
    pprint.pprint(par)

adata_file = par['adata']
output = par['output']
hvg = par['hvg']
scaling = par['scaling']

print('Read adata')
adata = sc.read(adata_file)

print('Prepare data')
if hvg > 0:
    # TODO: check that hvg value makes sense on dataset
    adata = hvg_batch(adata, batch_key='batch', target_genes=hvg, adataOut=True)

if scaling:
    adata = scale_batch(adata, batch='batch')

print('Integrate')
adata = runBBKNN(adata, batch='batch')

print('Save HDF5')
adata.uns['hvg'] = hvg
adata.uns['scaled'] = scaling
adata.write(output, compression='gzip')
