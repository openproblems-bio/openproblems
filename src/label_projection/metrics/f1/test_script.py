import os
from os import path
import subprocess

import scanpy as sc
import numpy as np

INPUT = "toy_baseline_pred_data.h5ad"
OUTPUT = "output.accuracy.h5ad"
AVG = "weighted"

print(">> Running f1 component")
out = subprocess.check_output([
    "./f1",
    "--input", INPUT,
    "--average", AVG,
    "--output", OUTPUT
]).decode("utf-8")

print(">> Checking whether file exists")
assert path.exists(OUTPUT)

print(">> Check that dataset fits expected API")
adata = sc.read_h5ad(OUTPUT)

# check id
assert "metric_id" in adata.uns
assert adata.uns["metric_id"] == "f1"
assert "metric_value" in adata.uns
assert type(adata.uns["metric_value"]) is np.float64

print(">> All tests passed successfully")
