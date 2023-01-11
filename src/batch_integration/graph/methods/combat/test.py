from os import path
import subprocess
import numpy as np
import scanpy as sc

np.random.seed(42)

method = 'combat'
output_file = method + '.h5ad'

print(">> Running script")
out = subprocess.check_output([
    "./" + method,
    "--input", 'processed.h5ad',
    "--hvg", 'False',
    "--scaling", 'False',
    "--output", output_file
]).decode("utf-8")

print(">> Checking whether file exists")
assert path.exists(output_file)

print('>> Checking API')
adata = sc.read(output_file, as_sparse='X')

assert 'dataset_id' in adata.uns
assert 'label' in adata.obs.columns
assert 'batch' in adata.obs.columns
assert 'highly_variable' in adata.var
assert 'counts' in adata.layers
assert 'logcounts' in adata.layers
assert 'logcounts_scaled' in adata.layers
assert 'X_pca' in adata.obsm
assert 'X_uni' in adata.obsm
assert 'uni' in adata.uns
assert 'uni_distances' in adata.obsp
assert 'uni_connectivities' in adata.obsp

assert 'connectivities' in adata.obsp
assert 'distances' in adata.obsp
assert 'hvg' in adata.uns
assert adata.uns['hvg'] == False
assert 'scaled' in adata.uns
assert adata.uns['scaled'] == False

unintegrated = sc.read('processed.h5ad', as_sparse='X')
assert len(unintegrated.X.data) != len(adata.X.data)
assert not np.any(np.not_equal(unintegrated.obsm['X_pca'], adata.obsm['X_pca']))

print(">> All tests passed successfully")
