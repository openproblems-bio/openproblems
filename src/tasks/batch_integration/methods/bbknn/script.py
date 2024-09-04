import sys
import anndata as ad
import scanpy as sc
import bbknn

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad',
    'annoy_n_trees': 10,
    'neighbors_within_batch': 3,
    'n_hvg': 2000,
}
meta = {
    'functionality_name': 'foo',
    'config': 'bar'
}
## VIASH END

sys.path.append(meta["resources_dir"])
from read_anndata_partial import read_anndata


print('Read input', flush=True)
adata = read_anndata(
    par['input'],
    X='layers/normalized',
    obs='obs',
    var='var',
    uns='uns'
)

if par['n_hvg']:
    print(f"Select top {par['n_hvg']} high variable genes", flush=True)
    idx = adata.var['hvg_score'].to_numpy().argsort()[::-1][:par['n_hvg']]
    adata = adata[:, idx].copy()
    sc.pp.pca(adata)

print('Run BBKNN', flush=True)
kwargs = dict(batch_key='batch', copy=True)
kwargs['annoy_n_trees'] = par['annoy_n_trees']
kwargs['neighbors_within_batch'] = par['neighbors_within_batch']

ad_bbknn = bbknn.bbknn(adata, **kwargs)

print("Store output", flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    obsp={
        'connectivities': ad_bbknn.obsp['connectivities'],
        'distances': ad_bbknn.obsp['distances'],
    },
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': meta['functionality_name'],
        'neighbors': ad_bbknn.uns['neighbors']
    }
)

print("Store outputs", flush=True)
output.write_h5ad(par['output'], compression='gzip')
