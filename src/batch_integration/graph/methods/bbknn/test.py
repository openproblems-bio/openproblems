from os import path
import subprocess
import numpy as np
import scanpy as sc

np.random.seed(42)

method = 'bbknn'
output_file = method + '.h5ad'

print(">> Running script")
out = subprocess.check_output([
    "./" + method,
    "--input", 'datasets_pancreas.h5ad',
    "--hvg", 'False',
    "--scaling", 'False',
    "--output", output_file
]).decode("utf-8")

print(">> Checking whether file exists")
assert path.exists(output_file)

print('>> Checking API')
adata = sc.read(output_file)

assert 'name' in adata.uns
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
assert not adata.uns['hvg']
assert 'scaled' in adata.uns
assert not adata.uns['scaled']

print(">> All tests passed successfully")
