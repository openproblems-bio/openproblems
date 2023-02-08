from os import path
import subprocess
import numpy as np
import anndata as ad


print(">> Running script", flush=True)

input_path = meta["resources_dir"] + "/pancreas/processed.h5ad"
output_path = "inegrated.h5ad"
cmd = [
    meta['executable'],
    "--input", input_path,
    "--output", output_path
]

print(">> Checking whether input file exists", flush=True)
assert path.exists(input_path)

print(">> Running script as test", flush=True)
subprocess.run(cmd, check=True)

print(">> Checking whether output file exists", flush=True)
assert path.exists(output_path)

print(">> Reading h5ad files", flush=True)
input = ad.read_h5ad(input_path)
output = ad.read_h5ad(output_path)

print(">> Checking whether predictions were added", flush=True)
assert 'dataset_id' in output.uns
assert 'X_pca' in output.obsm
assert 'X_emb' in output.obsm
assert 'normalization_id' in output.uns
assert 'method_id' in output.uns
assert meta['fuctionality_name'] == output.uns['method_id']

assert 'hvg' in output.uns
assert output.uns['hvg'] == False
assert 'scaled' in output.uns
assert output.uns['scaled'] == False

print(">> Checking whether data from input was copied properly to output", flush=True)
assert input.n_obs == output.n_obs
assert input.uns["dataset_id"] == output.uns["dataset_id"]

assert not np.any(np.not_equal(input.obsm['X_pca'], output.obsm['X_pca']))

print(">> All tests passed successfully")
