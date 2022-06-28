import subprocess
import scanpy as sc
from os import path

INPUT = "toy_preprocessed_data.h5ad"
OUTPUT = "output.scarchesallgenes.h5ad"

print(">> Running script as test")
out = subprocess.check_output([
    "./scarches_scanvi_all_genes",
    '--n_hidden', "32",
    '--n_layers', "1",
    '--n_latent', "10",
    '--max_epochs', "1",
    '--limit_train_batches', "10",
    '--limit_val_batches', "10",
    "--input", INPUT,
    "--output", OUTPUT
]).decode("utf-8")

print(">> Checking if output file exists")
assert path.exists(OUTPUT)

print(">> Checking if predictions were added")
adata = sc.read_h5ad(OUTPUT)
assert "celltype_pred" in adata.obs
assert "scarches_scanvi_all_genes" == adata.uns["method_id"]
