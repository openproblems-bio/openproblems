from os import path
import subprocess
import numpy as np
import anndata as ad
input_path = meta["resources_dir"] + "/pancreas/unintegrated.h5ad"
output_path = "embeddding.h5ad"
cmd = [
    meta['executable'],
    "--input", input_path,
    "--output", output_path
]
print(">> Checking whether input file exists", flush=True)
assert path.exists(input_path)

print(">> Running script as test", flush=True)
out = subprocess.run(cmd, stderr=subprocess.STDOUT).stdout
print(out)

print(">> Checking whether output file exists", flush=True)
assert path.exists(output_path)

print(">> Reading h5ad files", flush=True)
input = ad.read_h5ad(input_path)
output = ad.read_h5ad(output_path)
print(f"input: {input}", flush=True)
print(f"output: {output}", flush=True)

print(">> Checking whether predictions were added", flush=True)
assert 'dataset_id' in output.uns
assert 'X_pca' in output.obsm
assert 'X_emb' in output.obsm
assert 'normalization_id' in output.uns
assert 'method_id' in output.uns
assert meta['functionality_name'] == output.uns['method_id']
assert 'hvg' in output.uns

print(">> Checking whether data from input was copied properly to output", flush=True)
assert input.n_obs == output.n_obs
assert input.uns["dataset_id"] == output.uns["dataset_id"]
assert not np.any(np.not_equal(input.obsm['X_pca'], output.obsm['X_pca']))

print(">> All tests passed successfully")