## VIASH START
par = {
    'input': './src/batch_integration/graph/resources/mnn_pancreas.h5ad',
    'output': './src/batch_integration/graph/resources/ari_pancreas_mnn.tsv',
    'debug': True
}
## VIASH END

print('Importing libraries')
import pprint
import scanpy as sc
from scib.metrics.clustering import opt_louvain
from scib.metrics import ari

if par['debug']:
    pprint.pprint(par)

OUTPUT_TYPE = 'graph'
METRIC = 'ari'

adata_file = par['input']
output = par['output']

print('Read adata')
adata = sc.read(adata_file)
name = adata.uns['dataset_id']

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
score = ari(adata, group1='cluster', group2='label')

with open(output, 'w') as file:
    header = ['dataset', 'output_type', 'metric', 'value']
    entry = [name, OUTPUT_TYPE, METRIC, score]
    file.write('\t'.join(header) + '\n')
    file.write('\t'.join([str(x) for x in entry]))
