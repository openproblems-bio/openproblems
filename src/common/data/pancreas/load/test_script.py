import subprocess
import scanpy as sc
from os import path

URL = "https://ndownloader.figshare.com/files/24539828"
OUTPUT = "output.h5ad"

print(">> Running script")
out = subprocess.check_output([
    "./load_pancreas",
    "--url", URL,
    "--output", OUTPUT
]).decode("utf-8")

print(">> Checking whether file exists")
assert path.exists(OUTPUT)

print(">> Check that output fits expected API")
adata = sc.read_h5ad(OUTPUT)
# TODO: complete with API checks
assert "counts" not in adata.layers

print(">> All tests passed successfully")
