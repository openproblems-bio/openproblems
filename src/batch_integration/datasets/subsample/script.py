"""
This script subsets a downloaded dataset for testing all the batch integration modules
"""
## VIASH START
par = {
    'adata': 'src/common/dataset_loader/download/resources/pancreas.h5ad',
    'label': 'celltype',
    'batch': 'tech',
    'output': 'src/batch_integration/resources/data_loader_pancreas.h5ad',
    'debug': True
}
resources_dir = './src/batch_integration/datasets'
## VIASH END

print('Importing libraries')
import scanpy as sc
from pprint import pprint

if par['debug']:
    pprint(par)

adata_file = par['input']
label = par['label']
batch = par['batch']
output = par['output']
g2m_file = f'{resources_dir}/g2m_genes_tirosh_hm.txt'
s_file = f'{resources_dir}/s_genes_tirosh_hm.txt'

print(g2m_file)
import os
print(os.getcwd())
print(os.listdir())

print('Read adata')
adata = sc.read_h5ad(adata_file)

print('Get batch and label subsets')
head_batches = adata.obs[batch].unique().tolist()[0:3]
head_labels = adata.obs[label].unique().tolist()[0:2]

print('Get features subsets')
g2m_genes = [x.strip() for x in open(g2m_file).readlines()]
s_genes = [x.strip() for x in open(s_file).readlines()]

all_genes = adata.var.index.tolist()
cc_genes = [x for x in g2m_genes + s_genes if x in all_genes]
head_genes = list(set(cc_genes + all_genes[:100]))

print('Subset adata')
adata = adata[adata.obs[batch].isin(head_batches)]
adata = adata[adata.obs[label].isin(head_labels)]
adata = adata[:, head_genes]
sc.pp.subsample(adata, 0.3, random_state=42)

print(adata)

adata.write(output, compression='gzip')
