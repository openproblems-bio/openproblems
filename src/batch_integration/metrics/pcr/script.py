import anndata as ad
from scib.metrics import pcr_comparison

## VIASH START
par = {
    'input_integrated': 'resources_test/batch_integration/embedding/scvi.h5ad',
    'output': 'output.h5ad',
}

meta = {
    'functionality_name': 'foo',
}
## VIASH END

print('Read input', flush=True)
adata = ad.read_h5ad(par['input_integrated'])


print('preprocess data', flush=True)
adata.X = adata.layers['normalized']
adata_int = adata.copy()

print('compute score')
score = pcr_comparison(
    adata,
    adata_int,
    embed='X_emb',
    covariate='batch',
    verbose=False
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