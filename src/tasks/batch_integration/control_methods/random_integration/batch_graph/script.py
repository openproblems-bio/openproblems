import anndata as ad
import sys

## VIASH START

par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad'
}

meta = {
    'functionality_name': 'foo',
    'config': 'bar',
}

## VIASH END

# add helper scripts to path
sys.path.append(meta["resources_dir"])
from utils import _randomize_graph


print('Read input', flush=True)
adata = ad.read_h5ad(par['input'])

print('Randomize graph...', flush=True)
adata = _randomize_graph(
    adata,
    neighbors_key="knn",
    partition=adata.obs["batch"],
)

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['functionality_name']
adata.write_h5ad(par['output'], compression='gzip')
