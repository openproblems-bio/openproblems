import anndata as ad
import subprocess
from os import path

input_path = meta["resources_dir"] + "/pancreas/dataset.h5ad"
output_train_path = "output_train.h5ad"
output_test_path = "output_test.h5ad"

cmd = [
    meta['executable'],
    "--input", input_path,
    "--output_train", output_train_path,
    "--output_test", output_test_path,
]

print(">> Running script as test")
out = subprocess.run(cmd, check=True, capture_output=True, text=True)

print(">> Checking whether output file exists")
assert path.exists(output_train_path)
assert path.exists(output_test_path)

print(">> Reading h5ad files")
input = ad.read_h5ad(input_path)
output_train = ad.read_h5ad(output_train_path)
output_test = ad.read_h5ad(output_test_path)

print("input:", input)
print("output_train:", output_train)
print("output_test:", output_test)

print(">> Checking whether data from input was copied properly to output")
assert output_train.uns["dataset_id"] == input.uns["dataset_id"]
assert output_test.uns["dataset_id"] == input.uns["dataset_id"]

print(">> Check whether certain slots exist")
assert "counts" in output_train.layers
assert "counts" in output_test.layers

print(">> All checks succeeded!")