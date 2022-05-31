import subprocess
import scanpy as sc
from os import path


INPUT = "test_data.h5ad"
OUTPUT = "output.mv.h5ad"

print(">> Running script as test")
out = subprocess.check_output([
    "./majority_vote",
    "--input", INPUT,
    "--output", OUTPUT
]).decode("utf-8")

print(">> Checking if output file exists")
assert path.exists(OUTPUT)

print(">> Checking if predictions were added")
adata = sc.read_h5ad(OUTPUT)
assert "labels_pred" in adata.obs
assert "majority_vote" == adata.uns["method_id"]
