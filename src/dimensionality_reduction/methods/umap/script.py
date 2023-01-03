import anndata as ad
import scanpy as sc
import yaml

## VIASH START
par = {
    'input': 'resources_test/dimensionality_reduction/pancreas/train.h5ad',
    'output': 'reduced.h5ad',
    'n_pca': 50,
}
meta = {
    'functionality_name': 'umap',
    'config': 'src/dimensionality_reduction/methods/umap/config.vsh.yaml'
}
## VIASH END

print("Load input data")
input = ad.read_h5ad(par['input'])

print('Select top 1000 high variable genes')
n_genes = 1000
idx = input.var['hvg_score'].to_numpy().argsort()[::-1][:n_genes]
input = input[:, idx].copy()

print('Apply PCA with 50 dimensions')
input.obsm['X_pca_hvg'] = sc.tl.pca(input.layers['normalized'], n_comps=par['n_pca'], svd_solver="arpack")

print('Calculate a nearest-neighbour graph')
sc.pp.neighbors(input, use_rep="X_pca_hvg", n_pcs=par['n_pca'])

print("Run UMAP")
sc.tl.umap(input)
input.obsm["X_emb"] = input.obsm["X_umap"].copy()
del input.obsm["X_umap"]

print("Delete layers and var")
del input.layers
del input.var

print('Add method and normalization ID')
input.uns['method_id'] = meta['functionality_name']
with open(meta['config'], 'r') as config_file:
    config = yaml.safe_load(config_file)

input.uns['normalization_id'] = config['functionality']['info']['preferred_normalization']

print("Write output to file")
input.write_h5ad(par['output'], compression="gzip")