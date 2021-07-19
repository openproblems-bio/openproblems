## VIASH START
par = {
    'adata': '../resources/mnn.h5ad',
    'hvgs': 2000,
    'output': 'nmi.tsv',
    'debug': True
}
## VIASH END

print('Importing libraries')
import pprint
import scanpy as sc
from scIB.preprocessing import reduce_data
from scIB.clustering import opt_louvain
from scIB.metrics import nmi

if par['debug']:
    pprint.pprint(par)

METRIC = 'nmi'
EMBEDDING = 'X_pca'

adata_file = par['adata']
n_hvgs = par['hvgs']
n_hvgs = n_hvgs if n_hvgs > 0 else None
output = par['output']

print('Read adata')
adata = sc.read(adata_file)
name = adata.uns['name']

# preprocess adata object
print('preprocess adata')
reduce_data(
    adata,
    n_top_genes=n_hvgs,
    neighbors=True,
    use_rep=EMBEDDING,
    pca=True,
    umap=False
)

print('clustering')
opt_louvain(
    adata,
    label_key='label',
    cluster_key='cluster',
    plot=False,
    inplace=True,
    force=True
)

print('compute score')
score = nmi(adata, group1='cluster', group2='label')

with open(output, 'w') as file:
    header = ['dataset', 'output_type', 'hvg', 'metric', 'value']
    entry = [name, 'feature', n_hvgs, METRIC, score]
    file.write('\t'.join(header) + '\n')
    file.write('\t'.join([str(x) for x in entry]))
