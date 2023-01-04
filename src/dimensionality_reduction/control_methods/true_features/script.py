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

print("Load input data", flush=True)
input = ad.read_h5ad(par['input'])

print('Add method ID', flush=True)
input.uns['method_id'] = meta['functionality_name']

print('Create high dimensionally embedding with all features', flush=True)
input.obsm["X_emb"] = input.layers['counts'][:, :par['n_comps']].toarray()

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