import anndata as ad
import numpy as np

## VIASH START
par = {
    'input': 'resources_test/common/pancreas/dataset.h5ad',
    'output': 'output.h5ad',
}
meta = {
    'functionality_name': 'foo',
}
## VIASH END

print("Load input data")
adata = ad.read_h5ad(par['input'])

print('Create random embedding')
adata.obsm["X_emb"] = np.random.normal(0, 1, (adata.shape[0], 2))

# Update .uns
adata.uns['method_id'] = 'random_features'

print("Write output to file")
adata.write_h5ad(par['output'], compression="gzip")