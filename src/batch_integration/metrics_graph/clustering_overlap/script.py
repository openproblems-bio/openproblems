import anndata as ad
from scib.metrics.clustering import opt_louvain
from scib.metrics import ari, nmi

## VIASH START
par = {
    'input_integrated': 'resources_test/batch_integration/graph/bbknn.h5ad',
    'input_solution': 'resources_test/batch_integration/pancreas/solution.h5ad',
    'output': 'output.h5ad',
}

meta = {
    'functionality_name': 'foo'
}
## VIASH END

print('Read input', flush=True)
input = ad.read_h5ad(par['input_integrated'])
solution = ad.read_h5ad(par['input_solution'])

print('Transfer obs annotations', flush=True)
input.obs['label'] = solution.obs['label'][input.obs_names]

print('Run Louvain clustering', flush=True)
opt_louvain(
    input,
    label_key='label',
    cluster_key='cluster',
    plot=False,
    inplace=True,
    force=True
)

print('Compute ARI score', flush=True)
ari_score = ari(input, group1='cluster', group2='label')

print('Compute NMI score', flush=True)
nmi_score = nmi(input, group1='cluster', group2='label')

print("Create output AnnData object", flush=True)
output = ad.AnnData(
    uns={
        "dataset_id": input.uns['dataset_id'],
        'normalization_id': input.uns['normalization_id'],
        "method_id": input.uns['method_id'],
        "metric_ids": [ "ari", "nmi" ],
        "metric_values": [ ari_score, nmi_score ],
        "hvg": input.uns['hvg'],
        'output_type': input.uns['output_type'],
    }
)

if 'parent_method_id' in input.uns:
    output.uns['parent_method_id'] = input.uns['parent_method_id']

print("Write data to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")