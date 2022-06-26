import subprocess
import scanpy as sc
from os import path

INPUT = ["toy_normalized_log_cpm_data.h5ad", "toy_normalized_log_scran_pooling_data.h5ad"]
OUTPUT = ["output.lrlogcpm.h5ad", "output.lrscran.h5ad"]

for input, output in zip(INPUT, OUTPUT):
    print(">> Running script as test")
    out = subprocess.check_output([
        "./logistic_regression",
        "--input", input,
        "--max_iter", "100",
        "--output", output
    ]).decode("utf-8")

    print(">> Checking if output file exists")
    assert path.exists(output)

    print(">> Checking if predictions were added")
    adata = sc.read_h5ad(output)
    assert "celltype_pred" in adata.obs
    assert "logistic_regression" == adata.uns["method_id"]
