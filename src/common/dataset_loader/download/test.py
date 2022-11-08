from os import path
import subprocess
import scanpy as sc

name = "pancreas"
output = "dataset.h5ad"
url = "https://ndownloader.figshare.com/files/24539828"
obs_celltype = "celltype"
obs_batch = "tech"

print(">> Running script")
out = subprocess.check_output([
    meta["executable"],
    "--url", url,
    "--name", name,
    "--obs_celltype", obs_celltype,
    "--obs_batch", obs_batch,
    "--output", output
]).decode("utf-8")

print(">> Checking whether file exists")
assert path.exists(output)

print(">> Read output anndata")
adata = sc.read_h5ad(output)

print(">> Check that output fits expected API")
assert adata.X is not None
assert "counts" not in adata.layers
assert adata.uns["dataset_id"] == name
if obs_celltype:
    assert "celltype" in adata.obs.columns
if obs_batch:
    assert "batch" in adata.obs.columns

print(">> All tests passed successfully")
