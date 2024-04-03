import anndata as ad
import scanpy as sc
import scalex

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad',
    'n_hvg': 2000,
}
meta = {
    'functionality_name' : 'foo',
    'config': 'bar'
}
## VIASH END

print('Read input', flush=True)
adata = ad.read_h5ad(par['input'])

if par['n_hvg']:
    print(f"Select top {par['n_hvg']} high variable genes", flush=True)
    idx = adata.var['hvg_score'].to_numpy().argsort()[::-1][:par['n_hvg']]
    adata = adata[:, idx].copy()

print('Run SCALEX', flush=True)
adata.X = adata.layers['normalized']
adata = scalex.SCALEX(
    adata,
    batch_key="batch",
    ignore_umap=True,
    impute=adata.obs["batch"].cat.categories[0],
    processed=True,
    max_iteration=40,
    min_features=None,
    min_cells=None,
    n_top_features=0,
    outdir=None,
    gpu=0,
)

print("Store output", flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    layers={
        'corrected_counts': adata.layers["impute"],
    },
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': meta['functionality_name'],
    }
)

print("Write output to file", flush=True)
output.write_h5ad(par['output'], compression='gzip')
