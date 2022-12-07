import anndata as ad
import scanpy as sc
import yaml

## VIASH START
par = {
    'input': 'resources_test/common/pancreas/test.h5ad',
    'output': 'reduced.h5ad',
    'n_pca': 500,
}
meta = {
    'functionality_name': 'high_dim_pca',
}
## VIASH END

print("Load input data")
input = ad.read_h5ad(par['input'])

print('Add method and normalization ID')
input.uns['method_id'] = meta['functionality_name']
with open(meta['config'], 'r') as config_file:
    config = yaml.safe_load(config_file)

input.uns['normalization_id'] = config['functionality']['info']['preferred_normalization']

print('Create high dimensionally PCA embedding')
input.obsm["X_emb"] = sc.pp.pca(input.layers[input.uns['normalization_id']], n_comps=min(min(input.shape) - 1, par['n_pca']))

print("Write output to file")
input.write_h5ad(par['output'], compression="gzip")