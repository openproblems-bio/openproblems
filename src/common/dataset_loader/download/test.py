from os import path
import subprocess
import scanpy as sc

name = "pancreas"
output = "dataset.h5ad"
url = "https://ndownloader.figshare.com/files/24539828"
obs_celltype = "celltype"
obs_batch = "tech"

layer_counts_output = "foobar"

print(">> Running script")
out = subprocess.check_output([
    meta["executable"],
    "--url", url,
    "--name", name,
    "--obs_celltype", obs_celltype,
    "--obs_batch", obs_batch,
    "--layer_counts", "counts",
    "--layer_counts_output", layer_counts_output,
    "--output", output
]).decode("utf-8")

print(">> Checking whether file exists")
assert path.exists(output)

print(">> Read output anndata")
adata = sc.read_h5ad(output)

print(adata)

print(">> Check that output fits expected API")
if layer_counts_output is not None:
    assert adata.X is None
    assert layer_counts_output in adata.layers
else:
    assert adata.X is not None
    assert layer_counts_output not in adata.layers
assert adata.uns["dataset_id"] == name
if obs_celltype:
    assert "celltype" in adata.obs.columns
if obs_batch:
    assert "batch" in adata.obs.columns

print(">> All tests passed successfully")
