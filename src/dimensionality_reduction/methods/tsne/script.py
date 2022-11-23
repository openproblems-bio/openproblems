import anndata as ad
import scanpy as sc

## VIASH START
par = {
    'input': 'resources_test/common/pancreas/dataset.h5ad',
    'output': 'output.h5ad',
    'n_pca': 50,
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

print('Apply PCA with 50 dimensions')
adata.obsm['X_pca_hvg'] = sc.tl.pca(adata.layers['normalized'][:, idx], n_comps=par['n_pca'], svd_solver="arpack")

print('Run t-SNE')
sc.tl.tsne(adata, use_rep="X_pca_hvg", n_pcs=par['n_pca'])
adata.obsm["X_emb"] = adata.obsm["X_tsne"].copy()

# Update .uns
adata.uns['method_id'] = 'tsne'
adata.uns['normalization_id'] = 'log_cpm'
#del(adata.uns["pca_variance"])
#del(adata.uns["tsne"])
# Update .obsm/.varm 
#del(adata.obsm["X_tsne"])
#del(adata.varm["pca_loadings"])

print("Write output to file")
adata.write_h5ad(par['output'], compression="gzip")