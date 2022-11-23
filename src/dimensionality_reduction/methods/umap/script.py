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

print('Calculate a nearest-neighbour graph')
sc.pp.neighbors(adata, use_rep="X_pca_hvg", n_pcs=par['n_pca'])

print("Run UMAP")
sc.tl.umap(adata)

adata.obsm["X_emb"] = adata.obsm["X_umap"].copy()

# Update .uns
adata.uns['method_id'] = 'umap'
adata.uns['normalization_id'] = 'log_cpm'

print("Write output to file")
adata.write_h5ad(par['output'], compression="gzip")