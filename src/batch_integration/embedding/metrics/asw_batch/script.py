## VIASH START
par = {
    'adata': './src/batch_integration/embedding/resources/mnn_pancreas.h5ad',
    'output': './src/batch_integration/embedding/resources/asw_batch_pancreas_mnn.tsv',
    'debug': True
}
## VIASH END

print('Importing libraries')
import pprint
import scanpy as sc
from scib.metrics import silhouette_batch

if par['debug']:
    pprint.pprint(par)

OUTPUT_TYPE = 'embedding'
METRIC = 'asw_batch'
EMBEDDING = 'X_emb'

adata_file = par['adata']
output = par['output']

print('Read adata')
adata = sc.read(adata_file)
name = adata.uns['dataset_id']

print('compute score')
score = silhouette_batch(
    adata,
    batch_key='batch',
    group_key='label',
    embed=EMBEDDING,
)

with open(output, 'w') as file:
    header = ['dataset', 'output_type', 'metric', 'value']
    entry = [name, OUTPUT_TYPE, METRIC, score]
    file.write('\t'.join(header) + '\n')
    file.write('\t'.join([str(x) for x in entry]))
