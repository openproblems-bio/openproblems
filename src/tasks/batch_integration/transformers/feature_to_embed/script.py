import scanpy as sc
import yaml

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/integrated_feature.h5ad',
    'ouput': 'output.h5ad'
}
## VIASH END

print('Read input', flush=True)
adata= sc.read_h5ad(par['input'])


print('Run PCA', flush=True)
adata.obsm['X_emb'] = sc.pp.pca(
    adata.layers["corrected_counts"],
    n_comps=50,
    use_highly_variable=False,
    svd_solver='arpack',
    return_info=False
)

print('Store outputs', flush=True)
adata.write_h5ad(par['output'], compression='gzip')