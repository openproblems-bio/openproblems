import anndata as ad
import subprocess
from os import path

input_path = meta["resources_dir"] + "/pancreas/dataset.h5ad"
output_unintegrated_path = "unintegrated.h5ad"
cmd = [
    meta['executable'],
    "--input", input_path,
    "--output", output_unintegrated_path
]

print(">> Checking whether input file exists", flush=True)
assert path.exists(input_path)

print(">> Running script as test", flush=True)
out = subprocess.run(cmd, stderr=subprocess.STDOUT).stdout
print(out)

print(">> Checking whether output files exist", flush=True)
assert path.exists(output_unintegrated_path)

print(">> Reading h5ad files", flush=True)
input = ad.read_h5ad(input_path)
output_unintegrated = ad.read_h5ad(output_unintegrated_path)

print("input:", input, flush=True)
print("output_unintegrated:", output_unintegrated, flush=True)

print(">> Checking whether data from input was copied properly to output", flush=True)
assert input.n_obs == output_unintegrated.n_obs 
assert input.uns["dataset_id"] == output_unintegrated.uns["dataset_id"] 


print(">> Check whether certain slots exist", flush=True)
assert "counts" in output_unintegrated.layers
assert "normalized" in output_unintegrated.layers
assert 'hvg' in output_unintegrated.var

print("All checks succeeded!", flush=True)