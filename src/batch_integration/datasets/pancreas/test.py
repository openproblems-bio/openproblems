from os import path
import subprocess
import scanpy as sc

name = 'pancreas'
anndata_in = 'data_loader_pancreas.h5ad'
anndata_out = 'data_loader_pancreas.h5ad'

print('>> Running script')
out = subprocess.check_output([
    './pancreas',
    '--adata', anndata_in,
    '--label', 'celltype',
    '--batch', 'tech',
    '--hvgs', '2000',
    '--output', anndata_out
]).decode('utf-8')

print('>> Checking whether file exists')
assert path.exists(anndata_in)

print('>> Check that output fits expected API')
adata = sc.read_h5ad(anndata_in)
assert 'label' in adata.obs.columns
assert 'batch' in adata.obs.columns
assert 'logcounts' in adata.layers
assert 'X_pca' in adata.obsm
assert 'X_uni' in adata.obsm
assert 'uni' in adata.uns
assert 'uni_distances' in adata.obsp
assert 'uni_connectivities' in adata.obsp

print('>> All tests passed successfully')
