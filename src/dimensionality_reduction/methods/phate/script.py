import anndata as ad
from phate import PHATE
import scprep as sc
import yaml

## VIASH START
par = {
    'input': 'resources_test/common/pancreas/train.h5ad',
    'output': 'reduced.h5ad',
    'n_pca': 50,
    'g0': False,
    'log_cpm': False
}
meta = {
    'functionality_name': 'phate',
    'config': 'src/dimensionality_reduction/methods/phate/config.vsh.yaml'
}
## VIASH END

print("Load input data")
input = ad.read_h5ad(par['input'])

print("Run PHATE...")
gamma = 0 if par['g0'] else 1
print('... with gamma=' + str(gamma) + ' and...')
phate_op = PHATE(n_pca=par['n_pca'], verbose=False, n_jobs=-1, gamma=gamma)

with open(meta['config'], 'r') as config_file:
    config = yaml.safe_load(config_file)
input.uns['normalization_id'] = config['functionality']['info']['preferred_normalization']

if input.uns['normalization_id'] == 'sqrt_cpm':
    print('... using sqrt-CPM data')
    input.obsm["X_emb"] = phate_op.fit_transform(input.layers['normalized'])
elif input.uns['normalization_id'] == 'log_cpm':
    print('... using logCPM data')
    n_genes = 1000
    idx = input.var['hvg_score'].to_numpy().argsort()[::-1][:n_genes]
    input.obsm["X_emb"] = phate_op.fit_transform(input.layers['normalized'][:, idx])

print("Delete layers and var")
del input.layers
del input.var

print('Add method')
input.uns['method_id'] = meta['functionality_name']

print("Write output to file")
input.write_h5ad(par['output'], compression="gzip")