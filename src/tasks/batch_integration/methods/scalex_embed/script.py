import anndata as ad
import scanpy as sc
import scalex

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad',
    'hvg': True,
}
meta = {
    'functionality_name' : 'foo',
    'config': 'bar'
}
## VIASH END



print('Read input', flush=True)
adata = ad.read_h5ad(par['input'])

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
adata.obsm["X_emb"] = adata.obsm["latent"]

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['functionality_name']
adata.write_h5ad(par['output'], compression='gzip')
