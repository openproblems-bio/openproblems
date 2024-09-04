import anndata as ad
import mnnpy

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad',
    'n_hvg': 2000,
}
meta = {
    'functionality_name': 'foo',
    'config': 'bar'
}
## VIASH END

print('Read input', flush=True)
adata = ad.read_h5ad(par['input'])
adata.X = adata.layers['normalized']
del adata.layers['normalized']
del adata.layers['counts']

if par['n_hvg']:
    print(f"Select top {par['n_hvg']} high variable genes", flush=True)
    idx = adata.var['hvg_score'].to_numpy().argsort()[::-1][:par['n_hvg']]
    adata = adata[:, idx].copy()

print('Run mnn', flush=True)
split = []
batch_categories = adata.obs['batch'].cat.categories
for i in batch_categories:
    split.append(adata[adata.obs['batch'] == i].copy())
corrected, _, _ = mnnpy.mnn_correct(
        *split,
        batch_key='batch',
        batch_categories=batch_categories,
        index_unique=None
    )

print("Store outputs", flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': meta['functionality_name'],
    },
    layers={
        'corrected_counts': corrected.X,
    }
)


print("Store outputs", flush=True)
output.write_h5ad(par['output'], compression='gzip')
