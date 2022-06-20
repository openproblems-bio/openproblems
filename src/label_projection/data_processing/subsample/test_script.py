import subprocess
import scanpy as sc
from os import path

## VIASH START
## VIASH END

INPUT = f"{meta['resources_dir']}/pancreas/raw_data.h5ad"
OUTPUT = "toy_data.h5ad"

print(">> Runing script as test")
out = subprocess.check_output([
    "./" + meta["functionality_name"],
    "--input", INPUT,
    "--celltype_categories", "0:3",
    "--tech_categories", "0:-3:-2",
    "--output", OUTPUT
]).decode("utf-8")
print(">> Checking whether file exists")
assert path.exists(OUTPUT)

print(">> Check that test output fits expected API")
adata = sc.read_h5ad(OUTPUT)
assert (500, 443) == adata.X.shape, "processed result data shape {}".format(adata.X.shape)
