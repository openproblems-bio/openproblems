import anndata as ad
import scanpy as sc
from scib.metrics.clustering import cluster_optimal_resolution
from scib.metrics import ari, nmi

## VIASH START
par = {
    'input_integrated': 'resources_test/batch_integration/pancreas/integrated_graph.h5ad',
    'output': 'output.h5ad',
}

meta = {
    'functionality_name': 'foo'
}
## VIASH END

print('Read input', flush=True)
input_solution = ad.read_h5ad(par['input_solution'])
input_integrated = ad.read_h5ad(par['input_integrated'])

input_solution.obsp["connectivities"] = input_integrated.obsp["connectivities"]
input_solution.obsp["distances"] = input_integrated.obsp["distances"]

# TODO: if we don't copy neighbors over, the metric doesn't work
input_solution.uns["neighbors"] = input_integrated.uns["neighbors"]

print('Run optimal Leiden clustering', flush=True)
cluster_optimal_resolution(
    adata=input_solution,
    label_key='label',
    cluster_key='cluster',
    cluster_function=sc.tl.leiden,
)

print('Compute ARI score', flush=True)
ari_score = ari(input_solution, group1='cluster', group2='label')

print('Compute NMI score', flush=True)
nmi_score = nmi(input_solution, group1='cluster', group2='label')

print("Create output AnnData object", flush=True)
output = ad.AnnData(
    uns={
        "dataset_id": input_solution.uns['dataset_id'],
        'normalization_id': input_solution.uns['normalization_id'],
        "method_id": input_integrated.uns['method_id'],
        "metric_ids": [ "ari", "nmi" ],
        "metric_values": [ ari_score, nmi_score ]
    }
)

print("Write data to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")