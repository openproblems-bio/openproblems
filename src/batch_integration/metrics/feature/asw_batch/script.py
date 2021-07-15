## VIASH START
par = {
    'adata': '../resources/mnn.h5ad',
    'hvgs': 2000,
    'output': 'asw_batch.tsv',
    'debug': True
}
## VIASH END

print('Importing libraries')
import pprint
import scanpy as sc
from scIB.preprocessing import reduce_data
from scIB.metrics import silhouette_batch

if par['debug']:
    pprint.pprint(par)

METRIC = 'asw_batch'
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
    pca=True,
    neighbors=False,
    use_rep=EMBEDDING,
    umap=False
)

print('compute score')
_, sil_clus = silhouette_batch(
    adata,
    batch_key='batch',
    group_key='label',
    embed=EMBEDDING,
    verbose=False
)
score = sil_clus['silhouette_score'].mean()

with open(output, 'w') as file:
    header = ['dataset', 'output_type', 'hvg', 'metric', 'value']
    entry = [name, 'feature', n_hvgs, METRIC, score]
    file.write('\t'.join(header) + '\n')
    file.write('\t'.join([str(x) for x in entry]))
