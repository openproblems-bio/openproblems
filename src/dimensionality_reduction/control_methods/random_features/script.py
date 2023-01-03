import anndata as ad
import numpy as np
import yaml

## VIASH START
par = {
    'input': 'resources_test/common/pancreas/test.h5ad',
    'output': 'reduced.h5ad',
}
meta = {
    'functionality_name': 'random_features',
}
## VIASH END

print("Load input data", flush=True)
input = ad.read_h5ad(par['input'])

print('Add method and normalization ID', flush=True)
input.uns['method_id'] = meta['functionality_name']
with open(meta['config'], 'r') as config_file:
    config = yaml.safe_load(config_file)

input.uns['normalization_id'] = config['functionality']['info']['preferred_normalization']

print('Create random embedding', flush=True)
input.obsm["X_emb"] = np.random.normal(0, 1, (input.shape[0], 2))

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