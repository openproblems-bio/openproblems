import sys
import scanpy as sc
from scipy.sparse import csr_matrix

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad',
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


print('Run Combat', flush=True)
adata.X = sc.pp.combat(adata, key='batch', inplace=False)


print("Store output", flush=True)
output = sc.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': meta['functionality_name'],
    },
    layers={
        'corrected_counts': csr_matrix(adata.X),
    }
)

print("Store outputs", flush=True)
output.write_h5ad(par['output'], compression='gzip')
