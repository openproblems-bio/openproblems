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

print("Load input data", flush=True)
input = ad.read_h5ad(par['input'])

print('Select top 1000 high variable genes', flush=True)
n_genes = 1000
idx = input.var['hvg_score'].to_numpy().argsort()[::-1][:n_genes]
input = input[:, idx].copy()

print("Run UMAP...", flush=True)
if par['no_pca']:
    print('... using logCPM data', flush=True)
    input.obsm["X_emb"] = UMAP(densmap=True, random_state=42).fit_transform(input.layers['normalized'])
else:
    print('... after applying PCA with 50 dimensions to logCPM data', flush=True)
    input.obsm['X_pca_hvg'] = sc.tl.pca(input.layers['normalized'], n_comps=50, svd_solver="arpack")
    input.obsm["X_emb"] = UMAP(densmap=True, random_state=42).fit_transform(input.obsm['X_pca_hvg'])

print('Add method and normalization ID', flush=True)
input.uns['method_id'] = meta['functionality_name']
with open(meta['config'], 'r') as config_file:
    config = yaml.safe_load(config_file)

input.uns['normalization_id'] = config['functionality']['info']['preferred_normalization']

print('Copy data to new AnnData object', flush=True)
output = ad.AnnData(
    obs=input.obs[[]],
    uns={}
)
output.obsm['X_emb'] = input.obsm['X_emb']
output.uns['dataset_id'] = input.uns['dataset_id']
output.uns['normalization_id'] = input.uns['normalization_id']
output.uns['method_id'] = input.uns['method_id']

print("Write output to file", flush=True)
output.write_h5ad(par['output'], compression="gzip")