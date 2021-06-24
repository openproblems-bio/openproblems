from os import path
import subprocess
import scanpy as sc

name = "pbmc"
anndata_file = "pcmc.h5ad"

print(">> Running script")
out = subprocess.check_output([
    "./data_loader",
    "--url", "https://ndownloader.figshare.com/files/24974582",
    "--name", name,
    "--output", anndata_file
]).decode("utf-8")

print(">> Checking whether file exists")
assert path.exists(anndata_file)

print(">> Check that output fits expected API")
adata = sc.read_h5ad(anndata_file)
# TODO: complete with API checks
assert "counts" not in adata.layers
assert adata.uns["name"] == name

print(">> All tests passed successfully")
