import anndata as ad
from phate import PHATE
import scprep as sc
import yaml

## VIASH START
par = {
    'input': 'resources_test/dimensionality_reduction/pancreas/train.h5ad',
    'output': 'reduced.h5ad',
    'n_pca': 50,
    'g0': False,
}
meta = {
    'functionality_name': 'phate',
    'config': 'src/dimensionality_reduction/methods/phate/config.vsh.yaml'
}
## VIASH END

print("Load input data", flush=True)
input = ad.read_h5ad(par['input'])

print("Run PHATE...", flush=True)
gamma = 0 if par['g0'] else 1
print('... with gamma=' + str(gamma) + ' and...', flush=True)
phate_op = PHATE(n_pca=par['n_pca'], verbose=False, n_jobs=-1, gamma=gamma)

with open(meta['config'], 'r') as config_file:
    config = yaml.safe_load(config_file)
input.uns['normalization_id'] = config['functionality']['info']['preferred_normalization']

if input.uns['normalization_id'] == 'sqrt_cpm':
    print('... using sqrt-CPM data', flush=True)
    input.obsm["X_emb"] = phate_op.fit_transform(input.layers['normalized'])
elif input.uns['normalization_id'] == 'log_cpm':
    print('... using logCPM data', flush=True)
    n_genes = 1000
    idx = input.var['hvg_score'].to_numpy().argsort()[::-1][:n_genes]
    input = input[:, idx].copy()
    input.obsm["X_emb"] = phate_op.fit_transform(input.layers['normalized'])

print('Add method', flush=True)
input.uns['method_id'] = meta['functionality_name']

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