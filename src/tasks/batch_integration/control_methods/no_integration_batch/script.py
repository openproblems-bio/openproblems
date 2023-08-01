import scanpy as sc
import numpy as np
import yaml

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
input = sc.read_h5ad(par['input'])

print("process dataset", flush=True)
input.obsm["X_emb"] = np.zeros((input.shape[0], 50), dtype=float)
for batch in input.obs["batch"].unique():
    batch_idx = input.obs["batch"] == batch
    n_comps = min(50, np.sum(batch_idx))
    solver = "full" if n_comps == np.sum(batch_idx) else "arpack"
    # input.obsm["X_emb"][batch_idx, :n_comps] = sc.tl.pca(
    #     input[batch_idx],
    #     n_comps=n_comps,
    #     use_highly_variable=False,
    #     svd_solver=solver,
    #     copy=True,
    # ).obsm["X_pca"]

print("Store outputs", flush=True)
input.uns['method_id'] = meta['functionality_name']
input.write_h5ad(par['output'], compression='gzip')