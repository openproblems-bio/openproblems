import subprocess
import scanpy as sc
from os import path

INPUT = ["toy_normalized_log_cpm_data.h5ad", "toy_normalized_log_scran_pooling_data.h5ad"]
OUTPUT = ["result.logcpm.h5ad", "result.lrscran.h5ad"]
methods_params = {
        "./mlp": ["--hidden_layer_sizes", "20", "--max_iter", "100"],
        "./logistic_regression": ["--max_iter", "100"],
        "./knn_classifier": []
    }


_command = "./" + meta['functionality_name']
method_param = methods_params[_command]

for input, output in zip(INPUT, OUTPUT):
    default_params = ["--input", input, "--output", output]
    params = [_command] + default_params + method_param
    print(">> Running script as test")
    out = subprocess.check_output(params).decode("utf-8")

    print(">> Checking if output file exists")
    assert path.exists(output)

    print(">> Checking if predictions were added")
    adata = sc.read_h5ad(output)
    assert "celltype_pred" in adata.obs
    assert meta['functionality_name'] == adata.uns["method_id"]
