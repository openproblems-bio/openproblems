import anndata as ad
from scib.metrics import cell_cycle

## VIASH START
par = {
    'input_integrated': 'resources_test/batch_integration/pancreas/scvi.h5ad',
    'output': 'output.h5ad'
}

meta = {
    'functionality_name': 'foo'
}
## VIASH END

print('Read input', flush=True)
adata = ad.read_h5ad(par['input_integrated'])

adata.X = adata.layers['normalized']

adata_int = adata.copy()

translator = {
    "homo_sapiens": "human",
    "mus_musculus": "mouse"
}

print('compute score', flush=True)
score = cell_cycle(
    adata,
    adata_int,
    batch_key='batch',
    embed='X_emb',
    organism=translator[adata.uns['dataset_organism']]
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


print('Write data to file', flush=True)
output.write_h5ad(par['output'], compression='gzip')
