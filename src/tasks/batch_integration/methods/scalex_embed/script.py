import sys
import anndata as ad
import scalex

## VIASH START
par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad',
    'hvg': True,
}
meta = {
    'functionality_name' : 'foo',
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

print('Run SCALEX', flush=True)
adata = scalex.SCALEX(
    adata,
    batch_key="batch",
    ignore_umap=True,
    impute=adata.obs["batch"].cat.categories[0],
    processed=True,
    max_iteration=40,
    min_features=None,
    min_cells=None,
    n_top_features=0,
    outdir=None,
    gpu=0,
)
adata.obsm["X_emb"] = adata.obsm["latent"]

print("Store outputs", flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    layers={
        'corrected_counts': adata.layers["impute"],
    },
    obsm={
        'X_emb': adata.obsm['latent'],
    },
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'normalization_id': adata.uns['normalization_id'],
        'method_id': meta['functionality_name'],
    }
)

print("Write output to file", flush=True)
output.write_h5ad(par['output'], compression='gzip')
