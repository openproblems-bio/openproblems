import os
from os import path
import subprocess

import scanpy as sc

print(">> Running mnn")
out = subprocess.run([
    "./mnn",
    "--input", "sample_dataset.h5ad",
    "--output", "output.h5ad"
], stderr=subprocess.STDOUT)

if out.stdout:
    print(out.stdout)

if out.returncode:
    print(f"script: '{out.args}' exited with an error.")
    exit(out.returncode)

print(">> Checking whether file exists")
assert path.exists("output.h5ad")

print(">> Check that output fits expected API")
adata = sc.read_h5ad("output.h5ad")

assert "aligned" in adata.obsm
assert "mode2_aligned" in adata.obsm
assert adata.obsm["aligned"].shape[0] == adata.shape[0]
assert adata.obsm["mode2_aligned"].shape[0] == adata.obsm["mode2"].shape[0]
assert adata.obsm["aligned"].shape[1] == adata.obsm["mode2_aligned"].shape[1]

# check dataset id
assert "method_id" in adata.uns
assert adata.uns["method_id"] == "mnn"

print(">> All tests passed successfully")
