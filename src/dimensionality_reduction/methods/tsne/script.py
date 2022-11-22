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

adata = ad.read_h5ad(par['input'])

# t-SNE
sc.tl.tsne(adata, use_rep="X_pca", n_pcs=par['n_pca'])
adata.obsm["X_emb"] = adata.obsm["X_tsne"].copy()
del(adata.obsm["X_tsne"])
adata.uns["method_code_version"] = check_version("MulticoreTSNE")

adata.write_h5ad(par['output'], compression="gzip")
