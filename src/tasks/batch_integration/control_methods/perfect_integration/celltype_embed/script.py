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
from utils import _perfect_embedding
from read_anndata_partial import read_anndata


print('Read input', flush=True)
adata = read_anndata(
    par['input'],
    obs='obs',
    uns='uns'
)

print('Process data...', flush=True)
adata.obsm["X_emb"] = _perfect_embedding(partition=adata.obs["label"])

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['functionality_name']
adata.write_h5ad(par['output'], compression='gzip')
