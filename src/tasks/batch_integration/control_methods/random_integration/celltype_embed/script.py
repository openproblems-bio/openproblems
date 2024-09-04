import anndata as ad
import sys

## VIASH START

par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad',
}

meta = { 
    'functionality': 'foo',
    'config': 'bar'
}

## VIASH END
sys.path.append(meta["resources_dir"])
from utils import _randomize_features
from read_anndata_partial import read_anndata


print('Read input', flush=True)
adata = read_anndata(
    par['input'],
    obs='obs',
    obsm='obsm',
    uns='uns'
)

print('Process data...', flush=True)
adata.obsm["X_emb"] = _randomize_features(
    adata.obsm["X_pca"],
    partition=adata.obs["label"]
)

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['functionality_name']
adata.write_h5ad(par['output'], compression='gzip')