from os import path
import subprocess
import scanpy as sc

name = "pancreas"
output = "dataset.h5ad"
obs_celltype = "celltype"
obs_batch = "tech"

print(">> Running script", flush=True)
out = subprocess.run([
    meta["executable"],
    "--id", name,
    "--obs_celltype", obs_celltype,
    "--obs_batch", obs_batch,
    "--layer_counts", "counts",
    "--output", output],
    capture_output=True,
    text=True
)

print(out.stdout)

print(">> Checking whether file exists", flush=True)
assert path.exists(output)

print(">> Read output anndata", flush=True)
adata = sc.read_h5ad(output)

print(adata)

print(">> Check that output fits expected API", flush=True)
assert adata.X is None
assert "counts" in adata.layers
assert adata.uns["dataset_id"] == name
if obs_celltype:
    assert "celltype" in adata.obs.columns
if obs_batch:
    assert "batch" in adata.obs.columns

print(">> All tests passed successfully", flush=True)
