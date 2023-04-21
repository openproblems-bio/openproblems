import anndata as ad
import subprocess
from os import path

input_train_path = meta["resources_dir"] + "/pancreas/train.h5ad"
output_path = "output.h5ad"

cmd = [
    meta['executable'],
    "--input_train", input_train_path,
    "--output", output_path
]

print(">> Running script as test")
out = subprocess.run(cmd, check=True, capture_output=True, text=True)

print(">> Checking whether output file exists")
assert path.exists(output_path)

print(">> Reading h5ad files")
input_train = ad.read_h5ad(input_train_path)
output = ad.read_h5ad(output_path)
print("input_train:", input_train)
print("output:", output)

print(">> Checking whether predictions were added")
assert "denoised" in output.layers
assert meta['functionality_name'] == output.uns["method_id"]

print("Checking whether data from input was copied properly to output")
assert input_train.n_obs == output.n_obs
assert input_train.uns["dataset_id"] == output.uns["dataset_id"]

print("All checks succeeded!")
