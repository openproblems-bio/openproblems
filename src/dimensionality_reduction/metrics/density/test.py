import anndata as ad
import subprocess
from os import path

## VIASH START
meta = {
    'executable': './target/docker/dimensionality_reduction/density',
    'resources_dir': './resources_test/dimensionality_reduction/pancreas',
}
## VIASH END

input_reduced_path = meta["resources_dir"] + "/input/reduced.h5ad"
input_test_path = meta["resources_dir"] + "/input/test.h5ad"
output_path = "score.h5ad"
cmd = [
    meta['executable'],
    "--input_reduced", input_reduced_path,
    "--input_test", input_test_path,
    "--output", output_path,
]

print(">> Checking whether input files exist")
assert path.exists(input_reduced_path)
assert path.exists(input_test_path)

print(">> Running script as test")
out = subprocess.run(cmd, check=True, capture_output=True, text=True)

print(">> Checking whether output file exists")
assert path.exists(output_path)

print(">> Reading h5ad files")
input_reduced = ad.read_h5ad(input_reduced_path)
input_test = ad.read_h5ad(input_test_path)
output = ad.read_h5ad(output_path)

print("input reduced:", input_reduced)
print("input test:", input_test)
print("output:", output)

print(">> Checking whether metrics were added")
assert "metric_ids" in output.uns
assert "metric_values" in output.uns
assert meta['functionality_name'] in output.uns["metric_ids"]

print(">> Checking whether metrics are float")
assert isinstance(output.uns['metric_values'], float)

print(">> Checking whether data from input was copied properly to output")
assert input_reduced.uns["dataset_id"] == output.uns["dataset_id"]
assert input_reduced.uns["normalization_id"] == output.uns["normalization_id"]
assert input_reduced.uns["method_id"] == output.uns["method_id"]


print("All checks succeeded!")