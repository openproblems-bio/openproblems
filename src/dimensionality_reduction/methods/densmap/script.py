import anndata as ad
from umap import UMAP
import scanpy as sc

## VIASH START
par = {
    'input': 'resources_test/common/pancreas/dataset.h5ad',
    'output': 'output.h5ad',
    'no_pca': False,
}
meta = {
    'functionality_name': 'foo',
}
## VIASH END

print("Load input data")
adata = ad.read_h5ad(par['input'])

print('Select top 1000 high variable genes')
n_genes = 1000
idx = adata.var['hvg_score'].to_numpy().argsort()[::-1][:n_genes]

print("Run UMAP...")
if par['no_pca']:
    print('... using logCPM data')
    adata.obsm["X_emb"] = UMAP(densmap=True, random_state=42).fit_transform(adata.layers['normalized'][:, idx])
else:
    print('... after applying PCA with 50 dimensions to logCPM data')
    adata.obsm['X_pca_hvg'] = sc.tl.pca(adata.layers['normalized'][:, idx], n_comps=50, svd_solver="arpack")
    adata.obsm["X_emb"] = UMAP(densmap=True, random_state=42).fit_transform(adata.obsm['X_pca_hvg'])

print(adata.obsm['X_emb'][:10,:])

# Update .uns
adata.uns['method_id'] = 'densmap'
adata.uns['normalization_id'] = 'log_cpm'

print("Write output to file")
adata.write_h5ad(par['output'], compression="gzip")