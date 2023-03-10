import anndata as ad
from scib.metrics import cell_cycle

## VIASH START
par = {
    'input_integrated': 'resources_test/batch_integration/embedding/scvi.h5ad',
    'input_solution': 'resources_test/batch_integration/pancreas/solution.h5ad',
    'output': 'output.h5ad',
    'organism': 'human'
}

meta = {
    'functionality_name': 'foo'
}
## VIASH END

print('Read input', flush=True)
adata = ad.read_h5ad(par['input_integrated'])
adata_solution : ad.read_h5ad(par['input_solution'])


adata.X = adata.layers['normalized']

print('Transfer obs annotations', flush=True)
adata.obs['batch'] = adata_solution.obs['batch'][adata.obs_names]

adata_int = adata.copy()

print('compute score', flush=True)
score = cell_cycle(
    adata,
    adata_int,
    batch_key='batch',
    embed='X_emb',
    organism=par['organism']
)

print('Create output AnnData object', flush=True)
output = ad.AnnData(
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': adata.uns['method_id'],
        'metric_ids': [ meta['functionality_name'] ],
        'metric_values': [ score ],
        'hvg': adata.uns['hvg'],
        'output_type': adata.uns['output_type'],
    }
)

if 'parent_method_id' in adata.uns:
    output.uns['parent_method_id'] = adata.uns['parent_method_id']

print('Write data to file', flush=True)
output.write_h5ad(par['output'], compression='gzip')
