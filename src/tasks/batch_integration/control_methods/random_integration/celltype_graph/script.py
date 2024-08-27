import sys
import scanpy as sc

## VIASH START

par = {
    'input': 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
    'output': 'output.h5ad',
}

meta = { 
    'functionality': 'foo',
    'config': 'bar',
    "resources_dir": "src/tasks/batch_integration/control_methods/"
}

## VIASH END

# add helper scripts to path
sys.path.append(meta["resources_dir"])
from utils import _randomize_graph
from read_anndata_partial import read_anndata

print('Read input', flush=True)
adata = read_anndata(
    par['input'],
    obs='obs',
    obsp='obsp',
    uns='uns'
)

print("Process data...", flush=True)
adata = _randomize_graph(
    adata,
    neighbors_key="knn",
    partition=adata.obs["label"],
)

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['functionality_name']
adata.write_h5ad(par['output'], compression='gzip')
