from os import path
import subprocess
import scanpy as sc
import numpy as np

name = 'subsample'
anndata_in = 'pancreas.h5ad'
anndata_out = 'pancreas_sub.h5ad'

print('>> Running script')
n_hvgs = 100
out = subprocess.check_output([
    './subsample',
    '--input', anndata_in,
    '--label', 'celltype',
    '--batch', 'tech',
    '--output', anndata_out
]).decode('utf-8')

print('>> Checking whether file exists')
assert path.exists(anndata_out)

print('>> Check that output fits expected API')
adata = sc.read_h5ad(anndata_out)
assert 'dataset_id' in adata.uns

print('>> All tests passed successfully')
