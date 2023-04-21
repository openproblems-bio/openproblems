import anndata as ad
import subprocess
from os import path

input_denoised_path = meta["resources_dir"] + "/pancreas/magic.h5ad"
input_test_path = meta["resources_dir"] + "/pancreas/test.h5ad"
output_path = "output.h5ad"

cmd = [
    meta['executable'],
    "--input_denoised", input_denoised_path,
    "--input_test", input_test_path,
    "--output", output_path
]

print(">> Running script as test")
out = subprocess.run(cmd, check=True, capture_output=True, text=True)

print(">> Checking whether output file exists")
assert path.exists(output_path)

input_denoised = ad.read_h5ad(input_denoised_path)
input_test = ad.read_h5ad(input_test_path)
output = ad.read_h5ad(output_path)

print("Checking whether data from input was copied properly to output")
assert output.uns["dataset_id"] == input_denoised.uns["dataset_id"]
assert output.uns["method_id"] == input_denoised.uns["method_id"]
assert output.uns["metric_ids"] is not None
assert output.uns["metric_values"] is not None

# TODO: check whether the metric ids are all in .functionality.info

print("All checks succeeded!")