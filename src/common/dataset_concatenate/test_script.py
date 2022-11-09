import subprocess
import anndata as ad
from os import path

## VIASH START
## VIASH END

input_path = f"{meta['resources_dir']}/pancreas/dataset_cpm.h5ad"
output_path = "toy_data_concatenated.h5ad"

print(">> Runing script as test")
out = subprocess.check_output([
    meta["executable"],
    "--inputs", f"{input_path}:{input_path}",
    "--output", output_path
]).decode("utf-8")

print(">> Checking whether file exists")
assert path.exists(output_path)

print(">> Check that test output fits expected API")
input = ad.read_h5ad(input_path)
output = ad.read_h5ad(output_path)

assert output.n_obs == input.n_obs * 2
assert output.n_vars == input.n_vars
