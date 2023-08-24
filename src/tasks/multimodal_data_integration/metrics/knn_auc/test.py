import os
from os import path
import subprocess

import scanpy as sc
import numpy as np

print(">> Running knn_auc", flush=True)
out = subprocess.run([
    "./knn_auc",
    "--input", "sample_output.h5ad",
    "--output", "output.h5ad"
], stderr=subprocess.STDOUT)

if out.stdout:
    print(out.stdout)

if out.returncode:
    print(f"script: '{out.args}' exited with an error.")
    exit(out.returncode)

print(">> Checking whether file exists", flush=True)
assert path.exists("output.h5ad")

print(">> Check that dataset fits expected API", flush=True)
adata = sc.read_h5ad("output.h5ad")

# check id
assert "metric_id" in adata.uns
assert adata.uns["metric_id"] == "knn_auc"
assert "metric_value" in adata.uns
assert type(adata.uns["metric_value"]) is np.float64

print(">> All tests passed successfully", flush=True)
