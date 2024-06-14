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
from utils import _randomize_features


print('Read input', flush=True)
adata = sc.read_h5ad(par['input'])

print("Process data...", flush=True)
adata.layers['corrected_counts'] = _randomize_features(
    adata.layers["normalized"],
    partition=adata.obs["label"]
)

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['functionality_name']
adata.write_h5ad(par['output'], compression='gzip')
