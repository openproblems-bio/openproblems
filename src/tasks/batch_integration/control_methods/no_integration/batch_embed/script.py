import sys
import scanpy as sc
import numpy as np

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
from read_anndata_partial import read_anndata


print('Read input', flush=True)
adata = read_anndata(
    par['input'],
    X='layers/normalized',
    obs='obs',
    var='var',
    uns='uns'
)
adata.var["highly_variable"] = adata.var["hvg"]

print("Process dataset", flush=True)
adata.obsm["X_emb"] = np.zeros((adata.shape[0], 50), dtype=float)
for batch in adata.obs["batch"].unique():
    batch_idx = adata.obs["batch"] == batch
    n_comps = min(50, np.sum(batch_idx))
    solver = "full" if n_comps == np.sum(batch_idx) else "arpack"
    adata.obsm["X_emb"][batch_idx, :n_comps] = sc.tl.pca(
        adata[batch_idx].copy(),
        n_comps=n_comps,
        use_highly_variable=True,
        svd_solver=solver,
        copy=True,
    ).obsm["X_pca"]

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['functionality_name']
adata.write_h5ad(par['output'], compression='gzip')