import sys
import scanpy as sc

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/integrated_feature.h5ad',
    'ouput': 'output.h5ad'
}

meta = { 
    'functionality': 'foo',
    'config': 'bar'
}

## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata


print('Read input', flush=True)
adata = read_anndata(
    par['input'],
    X='layers/corrected_counts',
    obs='obs',
    var='var',
    uns='uns'
)


print('Run PCA', flush=True)
adata.obsm['X_emb'] = sc.pp.pca(
    adata.X,
    n_comps=50,
    use_highly_variable=False,  # Do we want to set this to True?
    svd_solver='arpack',
    return_info=False
)

print('Store outputs', flush=True)
adata.write_h5ad(par['output'], compression='gzip')