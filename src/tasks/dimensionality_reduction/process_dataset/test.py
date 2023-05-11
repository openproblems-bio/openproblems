import anndata as ad
import subprocess
from os import path

## VIASH START
meta = {
    'executable': './target/docker/dimensionality_reduction/',
    'resources_dir': './resources_test/common/',
}
## VIASH END

input_path = meta["resources_dir"] + "/input/dataset.h5ad"
output_train_path = "train.h5ad"
output_test_path = "test.h5ad"
cmd = [
    meta['executable'],
    "--input", input_path,
    "--output_train", output_train_path,
    "--output_test", output_test_path
]

print(">> Checking whether input file exists")
assert path.exists(input_path)

print(">> Running script as test")
out = subprocess.run(cmd, check=True, capture_output=True, text=True)

print(">> Checking whether output files exist")
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
assert input.n_obs == output_train.n_obs 
assert input.n_obs == output_test.n_obs
assert input.uns["dataset_id"] == output_train.uns["dataset_id"] 
assert input.uns["dataset_id"] == output_test.uns["dataset_id"]


print(">> Check whether certain slots exist")
assert "counts" in output_train.layers
assert "normalized" in output_train.layers
assert 'hvg_score' in output_train.var
assert "counts" in output_test.layers

print("All checks succeeded!")