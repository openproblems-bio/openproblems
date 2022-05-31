import os
from os import path
import subprocess

import scanpy as sc
import numpy as np

print(">> Running knn_auc")
out = subprocess.check_output([
    "./accuracy",
    "--input", "test_data.h5ad",
    "--output", "output.accuracy.h5ad"
]).decode("utf-8")

print(">> Checking whether file exists")
assert path.exists("output.accuracy.h5ad")

print(">> Check that dataset fits expected API")
adata = sc.read_h5ad("output.accuracy.h5ad")

# check id
assert "metric_id" in adata.uns
assert adata.uns["metric_id"] == "accuracy"
assert "metric_value" in adata.uns
assert type(adata.uns["metric_value"]) is np.float64

print(">> All tests passed successfully")
