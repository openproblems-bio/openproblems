from os import path
import subprocess
import scanpy as sc

name = "pancreas"
anndata_file = "dataset.h5ad"
url = "https://ndownloader.figshare.com/files/24974582"
obs_celltype = "celltype"
obs_batch = "tech"

print(">> Running script")
out = subprocess.check_output([
    f"./{meta['functionality_name']}",
    "--url", url,
    "--name", name,
    "--obs_celltype", obs_celltype,
    "--obs_batch", obs_batch,
    "--output", anndata_file
]).decode("utf-8")

print(">> Checking whether file exists")
assert path.exists(anndata_file)

print(">> Read output anndata")
adata = sc.read_h5ad(anndata_file)

print(">> Check that output fits expected API")
assert "counts" not in adata.layers
assert adata.uns["dataset_id"] == name
if obs_celltype:
    assert "celltype" in adata.obs
if obs_batch:
    assert "batch" in adata.obs

print(">> All tests passed successfully")
