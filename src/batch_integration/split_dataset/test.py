import os
import subprocess
import anndata as ad
import numpy as np

input_file = meta["resources_dir"] + '/pancreas/dataset.h5ad'
unintegrated_file = 'unintegrated.h5ad'
solution_file = 'solution.h5ad'
n_hvgs = 100

cmd_args = [
    meta["executable"],
    '--input', input_file,
    '--hvgs', str(n_hvgs),
    '--output_unintegrated', unintegrated_file,
    '--output_solution', solution_file
]
print('>> Running script')
subprocess.run(cmd_args, check=True)

print('>> Checking whether outputs exist')
assert os.path.exists(unintegrated_file)
assert os.path.exists(solution_file)

print('>> Read anndata files')
input = ad.read_h5ad(input_file)
unintegrated = ad.read_h5ad(unintegrated_file)
solution = ad.read_h5ad(solution_file)

print("input:", input)
print("unintegrated:", unintegrated)
print("solution:", solution)

print(">> Checking dimensions, make sure no cells were dropped")
assert input.n_obs == unintegrated.n_obs
assert input.n_obs == solution.n_obs
assert input.n_vars == unintegrated.n_vars
assert input.n_vars == solution.n_vars

print(">> Checking whether data from input was copied properly to output")
assert unintegrated.uns["dataset_id"] == input.uns["dataset_id"]
assert solution.uns["dataset_id"] == input.uns["dataset_id"]

print(">> Check output")
assert unintegrated.var['hvg'].dtype == 'bool'
assert unintegrated.var['hvg'].sum() == n_hvgs

print(">> Check whether certain slots exist")
# todo: use helper function for this

print('>> All tests passed successfully')
