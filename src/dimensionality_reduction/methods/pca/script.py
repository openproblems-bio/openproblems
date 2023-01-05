import anndata as ad
import scanpy as sc
import yaml

## VIASH START
par = {
    'input': 'resources_test/common/pancreas/train.h5ad',
    'output': 'reduced.h5ad',
}
meta = {
    'functionality_name': 'umap',
    'config': 'src/dimensionality_reduction/methods/PCA/config.vsh.yaml'
}
## VIASH END

print("Load input data", flush=True)
input = ad.read_h5ad(par['input'])

print('Select top 1000 high variable genes', flush=True)
n_genes = 1000
idx = input.var['hvg_score'].to_numpy().argsort()[::-1][:n_genes]

print('Apply PCA with 50 dimensions', flush=True)
input.obsm["X_emb"] = sc.tl.pca(input.layers['normalized'][:, idx], n_comps=50, svd_solver="arpack")[:, :2]

print('Add method ID', flush=True)
input.uns['method_id'] = meta['functionality_name']

print('Copy data to new AnnData object', flush=True)
output = ad.AnnData(
    obs=input.obs[[]],
    obsm={"X_emb": input.obsm["X_emb"]},
    uns={key: input.uns[key] for key in ["dataset_id", "normalization_id", "method_id"]}
)

print("Write output to file", flush=True)
output.write_h5ad(par['output'], compression="gzip")