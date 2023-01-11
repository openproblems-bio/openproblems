## VIASH START
par = {
    'adata': './src/batch_integration/embedding/resources/mnn_pancreas.h5ad',
    'output': './src/batch_integration/embedding/resources/cc_score_pancreas_mnn.tsv',
    'organism': 'human',
    'debug': True
}
## VIASH END

print('Importing libraries')
import pprint
import scanpy as sc
from scib.metrics import cell_cycle

if par['debug']:
    pprint.pprint(par)

OUTPUT_TYPE = 'embedding'
METRIC = 'cell_cycle_conservation'

adata_file = par['adata']
organism = par['organism']
output = par['output']

print('Read adata')
adata = sc.read(adata_file)
adata_int = adata.copy()
name = adata.uns['dataset_id']

print('compute score')
score = cell_cycle(
    adata,
    adata_int,
    batch_key='batch',
    embed='X_emb',
    organism=organism
)

with open(output, 'w') as file:
    header = ['dataset', 'output_type', 'metric', 'value']
    entry = [name, OUTPUT_TYPE, METRIC, score]
    file.write('\t'.join(header) + '\n')
    file.write('\t'.join([str(x) for x in entry]))
