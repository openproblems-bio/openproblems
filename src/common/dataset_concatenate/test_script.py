import subprocess
import scanpy as sc
from os import path

## VIASH START
## VIASH END
INPUT = f"{meta['resources_dir']}/pancreas/toy_data.h5ad"
INPUTS = f"{INPUT}:{INPUT}"
OUTPUT = "toy_data_concatenated.h5ad"

print(">> Runing script as test")
out = subprocess.check_output([
    "./" + meta["functionality_name"],
    "--inputs", INPUTS,
    "--output", OUTPUT
]).decode("utf-8")
print(">> Checking whether file exists")
assert path.exists(OUTPUT)

print(">> Check that test output fits expected API")
adata = sc.read_h5ad(OUTPUT)
assert (1000, 468) == adata.X.shape, "processed result data shape {}".format(adata.X.shape)
