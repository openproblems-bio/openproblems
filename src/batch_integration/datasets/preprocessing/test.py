import os
from os import path
import subprocess
import scanpy as sc
import numpy as np

name = 'preprocessing'
anndata_in = meta["resources_dir"] + '/pancreas/dataset.h5ad'
anndata_out = 'datasets_pancreas.h5ad'

print('>> Running script')
n_hvgs = 100
out = subprocess.check_output([
    './preprocessing',
    '--input', anndata_in,
    '--label', 'celltype',
    '--batch', 'tech',
    '--hvgs', str(n_hvgs),
    '--output', anndata_out
]).decode('utf-8')

print('>> Checking whether file exists')
assert path.exists(anndata_out)

print('>> Check that output fits expected API')
adata = sc.read_h5ad(anndata_out)
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

assert adata.var['highly_variable'].dtype == 'bool'
assert adata.var['highly_variable'].sum() == n_hvgs
assert -0.0000001 <= np.mean(adata.layers['logcounts_scaled']) <= 0.0000001
assert 0.75 <= np.var(adata.layers['logcounts_scaled']) <= 1

print('>> All tests passed successfully')
