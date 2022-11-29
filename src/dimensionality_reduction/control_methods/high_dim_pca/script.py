import anndata as ad
import scanpy as sc

## VIASH START
par = {
    'input': 'resources_test/common/pancreas/dataset.h5ad',
    'output': 'output.h5ad',
}
meta = {
    'functionality_name': 'foo',
}
## VIASH END

print("Load input data")
adata = ad.read_h5ad(par['input'])

print('Create high dimensionally PCA embedding')
sc.pp.pca(adata.layers['counts'], n_comps=min(min(adata.shape), par['n_pca']))
adata.obsm["X_emb"] = adata.obsm["X_pca"]

# Update .uns
adata.uns['method_id'] = 'high_dim_pca'

print("Write output to file")
adata.write_h5ad(par['output'], compression="gzip")