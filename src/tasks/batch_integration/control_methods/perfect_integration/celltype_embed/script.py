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


print('Read input', flush=True)
adata = ad.read_h5ad(par['input'])

print('Process data...', flush=True)
adata.obsm["X_emb"] = _perfect_embedding(partition=adata.obs["label"])

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['functionality_name']
adata.write_h5ad(par['output'], compression='gzip')
