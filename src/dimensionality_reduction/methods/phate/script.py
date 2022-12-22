import anndata as ad
from phate import PHATE

## VIASH START
par = {
    'input': 'resources_test/dimensionality_reduction/pancreas/train.h5ad',
    'output': 'reduced.h5ad',
    'n_pca': 50,
    'gamma': 1,
    'num_hvg_genes': None
}
meta = {
    'functionality_name': 'phate',
    'config': 'src/dimensionality_reduction/methods/phate/config.vsh.yaml'
}
## VIASH END

print("Load input data", flush=True)
input = ad.read_h5ad(par['input'])

print("Run PHATE", flush=True)
phate_op = PHATE(n_pca=par['n_pca'], verbose=False, n_jobs=-1, gamma=par['gamma'])
X_mat = input.layers['normalized']

if par["num_hvg_genes"]:
    print("Subsetting to hvg genes", flush=True)
    num_features = par["num_hvg_genes"]
    hvg_idx = input.var['hvg_score'].to_numpy().argsort()[::-1][:num_features]
    X_mat = X_mat[:, hvg_idx]

# store embedding
input.obsm["X_emb"] = phate_op.fit_transform(X_mat)

print("Delete layers and var", flush=True)
del input.layers
del input.var

print('Add method', flush=True)
input.uns['method_id'] = meta['functionality_name']

print("Write output to file", flush=True)
input.write_h5ad(par['output'], compression="gzip")