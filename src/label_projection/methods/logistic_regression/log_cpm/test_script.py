import subprocess
import scanpy as sc
from os import path

INPUT = "toy_preprocessed_data.h5ad"
OUTPUT = "output.lrlogcpm.h5ad"

print(">> Running script as test")
out = subprocess.check_output([
    "./logistic_regression_log_cpm",
    "--input", INPUT,
    "--max_iter", "100",
    "--output", OUTPUT
]).decode("utf-8")

print(">> Checking if output file exists")
assert path.exists(OUTPUT)

print(">> Checking if predictions were added")
adata = sc.read_h5ad(OUTPUT)
assert "labels_pred" in adata.obs
assert "logistic_regression_log_cpm" == adata.uns["method_id"]
