import subprocess
import scanpy as sc
from os import path

INPUT = "toy_preprocessed_data.h5ad"
OUTPUT = "output.scanviallgenes.h5ad"

print(">> Running script as test")
out = subprocess.check_output([
    "./scanvi_hvg",
    '--n_hidden', "32",
    '--n_layers', "1",
    '--n_latent', "10",
    '--n_top_genes', "2000",
    '--spane', "0.8",
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
assert "labels_pred" in adata.obs
assert "scanvi_hvg" == adata.uns["method_id"]
