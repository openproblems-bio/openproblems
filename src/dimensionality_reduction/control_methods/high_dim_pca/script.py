import anndata as ad
import scanpy as sc

## VIASH START
par = {
    'input': 'resources_test/common/pancreas/dataset.h5ad',
    'output': 'output.h5ad',
    'n_pca': 500,
}
meta = {
    'functionality_name': 'foo',
}
## VIASH END

print("Load input data")
adata = ad.read_h5ad(par['input'])

print('Create high dimensionally PCA embedding')
adata.obsm["X_emb"] = sc.pp.pca(adata.layers['counts'], n_comps=min(min(adata.shape) - 1, par['n_pca']))

# Update .uns
adata.uns['method_id'] = 'high_dim_pca'

print("Write output to file")
adata.write_h5ad(par['output'], compression="gzip")