import anndata as ad
from phate import PHATE
import scprep as sc
# import yaml

## VIASH START
par = {
    'input': 'resources_test/common/pancreas/dataset.h5ad',
    'output': 'output.h5ad',
    'n_pca': 50,
    'g0': False,
    'hvg': False
}
meta = {
    'functionality_name': 'foo',
    'config': 'src/dimensionality_reduction/methods/phate/config.vsh.yaml'
}
## VIASH END
# with open(meta['config'], 'r') as config_file:
#     config = yaml.safe_load(config_file)

# config['functionality']['info']['preferred_normalization']
# print(meta)

print("Load input data")
adata = ad.read_h5ad(par['input'])

print("Run PHATE...")
gamma = 0 if par['g0'] else 1
print('... with gamma=' + str(gamma) + ' and...')
phate_op = PHATE(n_pca=par['n_pca'], verbose=False, n_jobs=-1, gamma=gamma)

if par['hvg']:
    print('... using logCPM data')
    n_genes = 1000
    idx = adata.var['hvg_score'].to_numpy().argsort()[::-1][:n_genes]
    adata.obsm["X_emb"] = phate_op.fit_transform(adata.layers['normalized'][:, idx])
    adata.uns['normalization_id'] = 'log_cpm'
else:
    print('... using sqrt-CPM data')
    adata.layers['sqrt_cpm'] = sc.transform.sqrt(adata.layers['normalized'].expm1())
    adata.obsm["X_emb"] = phate_op.fit_transform(adata.layers['sqrt_cpm'])
    adata.uns['normalization_id'] = 'sqrt_cpm'

# Update .uns
adata.uns['method_id'] = 'phate'

print("Write output to file")
adata.write_h5ad(par['output'], compression="gzip")