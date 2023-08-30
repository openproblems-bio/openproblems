import anndata as ad
from scib.metrics import isolated_labels_f1

## VIASH START
par = {
    'input_integrated': 'resources_test/batch_integration/pancreas/integrated_embedding.h5ad',
    'output': 'output.h5ad',
}

meta = {
    'functionality_name': 'foo',
}
## VIASH END

print('Read input', flush=True)
input_solution = ad.read_h5ad(par['input_solution'])
input_integrated = ad.read_h5ad(par['input_integrated'])

input_solution.obsp["connectivities"] = input_integrated.obsp["connectivities"]
input_solution.obsp["distances"] = input_integrated.obsp["distances"]

# TODO: if we don't copy neighbors over, the metric doesn't work
input_solution.uns["neighbors"] = input_integrated.uns["neighbors"]

print('compute score', flush=True)
score = isolated_labels_f1(
    input_solution,
    label_key='label',
    batch_key='batch',
    embed=None,
    iso_threshold=None,
    verbose=True,
)
print(score, flush=True)

print('Create output AnnData object', flush=True)
output = ad.AnnData(
    uns={
        'dataset_id': input_solution.uns['dataset_id'],
        'normalization_id': input_solution.uns['normalization_id'],
        'method_id': input_integrated.uns['method_id'],
        'metric_ids': [ meta['functionality_name'] ],
        'metric_values': [ score ]
    }
)

print('Write data to file', flush=True)
output.write_h5ad(par['output'], compression='gzip')