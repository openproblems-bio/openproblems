import subprocess
import scanpy as sc
from os import path

INPUT = "test_data.h5ad"
OUTPUT = "preprocessed.h5ad"
METHODS = ["batch", "random", "random_with_noise"]

print(">> Runing script as test")
out = subprocess.check_output([
    "./preprocess_pancreas",
    "--input", INPUT,
    "--test", "True",
    "--output", OUTPUT
]).decode("utf-8")
print(">> Checking whether file exists")
assert path.exists(OUTPUT)

print(">> Check that test output fits expected API")
adata = sc.read_h5ad(OUTPUT)
assert (500, 441) == adata.X.shape, "processed result data shape {}".format(adata.X.shape)

for method in METHODS:
    print(">> Running script for {} method".format(method))
    out = subprocess.check_output([
        "./preprocess_pancreas",
        "--input", INPUT,
        "--method", method,
        "--output", OUTPUT
    ]).decode("utf-8")
    print(">> Checking whether file exists")
    assert path.exists(OUTPUT)
    print(">> Check that test output fits expected API")
    adata = sc.read_h5ad(OUTPUT)
    assert (16382, 18771) == adata.X.shape, "processed result data shape {}".format(adata.X.shape)
    assert "batch" in adata.obs
    assert "is_train" in adata.obs
