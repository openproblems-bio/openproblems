import subprocess
import scanpy as sc
from os import path

INPUT = "toy_preprocessed_data.h5ad"
OUTPUT = "result.scvimethod.h5ad"

methods_params = {
        "./scarches_scanvi_all_genes": ['--n_hidden', "32", '--n_layers', "1",
                                        '--n_latent', "10", '--max_epochs', "1",
                                        '--limit_train_batches', "10",
                                        '--limit_val_batches', "10"],
        "./scarches_scanvi_hvg": ['--n_hidden', "32", '--n_layers', "1",
                                  '--n_latent', "10",  '--span', "0.8",
                                  '--n_top_genes', "2000", '--max_epochs', "1",
                                  '--limit_train_batches', "10", '--limit_val_batches', "10"],
        "./scanvi_hvg": ['--n_hidden', "32", '--n_layers', "1",
                         '--n_latent', "10",  '--span', "0.8",
                         '--n_top_genes', "2000", '--max_epochs', "1",
                         '--limit_train_batches', "10", '--limit_val_batches', "10"],
        "./scanvi_all_genes": ['--n_hidden', "32", '--n_layers', "1",
                               '--n_latent', "10", '--max_epochs', "1",
                               '--limit_train_batches', "10", '--limit_val_batches', "10"]
    }


_command = "./" + meta['functionality_name']
method_param = methods_params[_command]

default_params = ["--input", INPUT, "--output", OUTPUT]
params = [_command] + default_params + method_param

print(">> Running script as test")
out = subprocess.check_output(params).decode("utf-8")

print(">> Checking if output file exists")
assert path.exists(OUTPUT)

print(">> Checking if predictions were added")
adata = sc.read_h5ad(OUTPUT)
assert "celltype_pred" in adata.obs
assert meta['functionality_name'] == adata.uns["method_id"]
