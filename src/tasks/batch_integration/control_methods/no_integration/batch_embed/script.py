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

print('Read input', flush=True)
adata = sc.read_h5ad(par['input'])
adata.X = adata.layers["normalized"]
adata.var["highly_variable"] = adata.var["hvg"]

print("Process dataset", flush=True)
adata.obsm["X_emb"] = np.zeros((adata.shape[0], 50), dtype=float)
for batch in adata.obs["batch"].unique():
    batch_idx = adata.obs["batch"] == batch
    n_comps = min(50, np.sum(batch_idx))
    solver = "full" if n_comps == np.sum(batch_idx) else "arpack"
    adata.obsm["X_emb"][batch_idx, :n_comps] = sc.tl.pca(
        adata[batch_idx],
        n_comps=n_comps,
        use_highly_variable=True,
        svd_solver=solver,
        copy=True,
    ).obsm["X_pca"]

print("Store outputs", flush=True)
adata.uns['method_id'] = meta['functionality_name']
adata.write_h5ad(par['output'], compression='gzip')