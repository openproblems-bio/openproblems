import anndata as ad
from umap import UMAP
import scanpy as sc
import yaml

## VIASH START
par = {
    'input': 'resources_test/common/pancreas/train.h5ad',
    'output': 'reduced.h5ad',
    'no_pca': False,
}
meta = {
    'functionality_name': 'densmap',
    'config': 'src/dimensionality_reduction/methods/densmap/config.vsh.yaml'
}
## VIASH END

print("Load input data")
input = ad.read_h5ad(par['input'])

print('Select top 1000 high variable genes')
n_genes = 1000
idx = input.var['hvg_score'].to_numpy().argsort()[::-1][:n_genes]

print("Run UMAP...")
if par['no_pca']:
    print('... using logCPM data')
    input.obsm["X_emb"] = UMAP(densmap=True, random_state=42).fit_transform(input.layers['normalized'][:, idx])
else:
    print('... after applying PCA with 50 dimensions to logCPM data')
    input.obsm['X_pca_hvg'] = sc.tl.pca(input.layers['normalized'][:, idx], n_comps=50, svd_solver="arpack")
    input.obsm["X_emb"] = UMAP(densmap=True, random_state=42).fit_transform(input.obsm['X_pca_hvg'])

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