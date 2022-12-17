import anndata as ad
import subprocess
from os import path

## VIASH START
meta = {
    'executable': './target/docker/dimensionality_reduction/umap',
    'resources_dir': './resources_test/common/',
}
## VIASH END

input_path = meta["resources_dir"] + "/input/train.h5ad"
output_path = "reduced.h5ad"
cmd = [
    meta['executable'],
    "--input", input_path,
    "--output", output_path
]

print(">> Checking whether input file exists")
assert path.exists(input_path)

print(">> Running script as test")
out = subprocess.run(cmd, check=True, capture_output=True, text=True)

print(">> Checking whether output file exists")
assert path.exists(output_path)

print(">> Reading h5ad files")
input = ad.read_h5ad(input_path)
output = ad.read_h5ad(output_path)

print("input:", input)
print("output:", output)

print(">> Checking whether predictions were added")
assert "X_emb" in output.obsm
assert meta['functionality_name'] == output.uns["method_id"]
assert 'normalization_id' in output.uns

print(">> Checking whether data from input was copied properly to output")
assert input.n_obs == output.n_obs
assert input.uns["dataset_id"] == output.uns["dataset_id"]

print("All checks succeeded!")