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

if par["num_hvg_genes"] and input.uns['normalization_id'] == 'log_cpm':
    print("Subsetting to hvg genes", flush=True)
    num_features = par["num_hvg_genes"]
    hvg_idx = input.var['hvg_score'].to_numpy().argsort()[::-1][:num_features]
    X_mat = X_mat[:, hvg_idx]

# store embedding
input.obsm["X_emb"] = phate_op.fit_transform(X_mat)

print('Add method', flush=True)
input.uns['method_id'] = meta['functionality_name']

print('Copy data to new AnnData object', flush=True)
output = ad.AnnData(
    obs=input.obs[[]],
    obsm={"X_emb": input.obsm["X_emb"]},
    uns={key: input.uns[key] for key in ["dataset_id", "normalization_id", "method_id"]}
)

print("Write output to file", flush=True)
output.write_h5ad(par['output'], compression="gzip")