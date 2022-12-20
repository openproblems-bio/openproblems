import anndata as ad
import scanpy as sc
import yaml

## VIASH START
par = {
    'input': 'resources_test/common/pancreas/test.h5ad',
    'output': 'reduced.h5ad',
    'n_comps': 100,
}
meta = {
    'functionality_name': 'true_features',
}
## VIASH END

print("Load input data")
input = ad.read_h5ad(par['input'])

print('Add method and normalization ID')
input.uns['method_id'] = meta['functionality_name']
with open(meta['config'], 'r') as config_file:
    config = yaml.safe_load(config_file)

input.uns['normalization_id'] = config['functionality']['info']['preferred_normalization']

print('Create high dimensionally embedding with all features')
input.obsm["X_emb"] = input.layers['counts'][:, :par['n_comps']].toarray()

print("Write output to file")
input.write_h5ad(par['output'], compression="gzip")