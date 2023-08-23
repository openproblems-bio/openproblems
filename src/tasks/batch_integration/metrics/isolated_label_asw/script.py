import anndata as ad
from scib.metrics import isolated_labels_asw

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
adata = ad.read_h5ad(par['input_integrated'])

print('preprocess data', flush=True)
adata.X = adata.layers['normalized']
adata_int = adata.copy()

print('compute score')

score = isolated_labels_asw(
    adata,
    label_key='label',
    batch_key='batch',
    embed='X_emb',
    iso_threshold=None,
    verbose=True,
)
print(score, flush=True)

print('Create output AnnData object', flush=True)
output = ad.AnnData(
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': adata.uns['method_id'],
        'metric_ids': [ meta['functionality_name'] ],
        'metric_values': [ score ]
    }
)

print('Write data to file', flush=True)
output.write_h5ad(par['output'], compression='gzip')