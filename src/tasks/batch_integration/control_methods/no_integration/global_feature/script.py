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

# no processing, subset matrix to highly variable genes
adata_hvg = adata[:, adata.var["hvg"]].copy()
adata.layers['corrected_counts'] = adata_hvg.X.copy()

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['functionality_name']
adata.write_h5ad(par['output'], compression='gzip')
