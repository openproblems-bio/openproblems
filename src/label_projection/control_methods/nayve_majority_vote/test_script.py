import subprocess
import scanpy as sc
from os import path


INPUT = "toy_preprocessed_data.h5ad"
OUTPUT = "output.mv.h5ad"

print(">> Running script as test")
out = subprocess.check_output([
    "./" + meta["functionality_name"],
    "--input", INPUT,
    "--output", OUTPUT
]).decode("utf-8")

print(">> Checking if output file exists")
assert path.exists(OUTPUT)

print(">> Checking if predictions were added")
adata = sc.read_h5ad(OUTPUT)
assert "celltype_pred" in adata.obs
assert meta["functionality_name"] == adata.uns["method_id"]
