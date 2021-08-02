from os import path
import subprocess
import numpy as np
import scanpy as sc

np.random.seed(42)

method = 'bbknn'
output_file = method + '_scaled_hvg.h5ad'

print(">> Running script")
out = subprocess.check_output([
    "./" + method,
    "--adata", 'datasets_pancreas.h5ad',
    "--hvg", '100',
    "--scaling", 'True',
    "--output", output_file
]).decode("utf-8")

print(">> Checking whether file exists")
assert path.exists(output_file)

print('>> Checking API')
adata = sc.read(output_file)
assert 'connectivities' in adata.obsp
assert 'distances' in adata.obsp
assert 'hvg' in adata.uns
assert adata.uns['hvg'] == 100
assert 'scaled' in adata.uns
print(isinstance(adata.uns['scaled'], bool))
print(type(adata.uns['scaled']).__name__())
assert adata.uns['scaled'] == True

print(">> All tests passed successfully")
