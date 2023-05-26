import subprocess
from os import path
import json

input_path = meta["resources_dir"] + "/pancreas/dataset.h5ad"
input_correct_schema = "anndata_correct.yaml"
input_error_schema = "anndata_error.yaml"
output_checks = "checks.json"
output_path = "output.h5ad"
output_error_checks = "error_checks.json"
output_error_path = "error_output.h5ad"



with open(input_correct_schema, "w") as f:
    f.write('''
type: file
description: "A preprocessed dataset"
example: "preprocessed.h5ad"
info:
  short_description: "Preprocessed dataset"
  slots:
    layers: 
      - type: integer
        name: counts
        description: Raw counts
    uns:
      - type: string
        name: dataset_id
        description: "A unique identifier for the dataset"
    ''')


with open(input_error_schema, "w") as f:
    f.write('''
type: file
description: "A preprocessed dataset"
example: "preprocessed.h5ad"
info:
  short_description: "Preprocessed dataset"
  slots:
    layers: 
      - type: integer
        name: counts
        description: Raw counts
    uns:
      - type: string
        name: dataset_id
        description: "A unique identifier for the dataset"
      - type: string
        name: error_test
        description: "A made up uns variable to test if error is picked up"
    ''')


cmd = [
    meta['executable'],
    "--input", input_path,
    "--schema", input_correct_schema,
    "--checks", output_checks,
    "--output", output_path,
]

print(">> Running script as test", flush=True)
out = subprocess.run(cmd, stderr=subprocess.STDOUT)

if out.stdout:
    print(out.stdout)

if out.returncode:
    print(f"script: '{cmd}' exited with an error.")
    exit(out.returncode)

print(">> Checking whether output file exists", flush=True)
assert path.exists(output_checks), "Output checks file does not exist"
assert path.exists(output_path), "Output path does not exist"

print(">> Reading json file", flush=True)
with open(output_checks, 'r') as f:
    out = json.load(f)
    print(out)


# Check if an incomplete h5ad is captured
cmd_error = [
    meta['executable'],
    "--input", input_path,
    "--schema", input_error_schema,
    "--stop_on_error", 'true',
    "--checks", output_error_checks,
    "--output", output_error_path,
]

print(">> Running script as test", flush=True)
out_error = subprocess.run(cmd_error)

print(">> Checking whether output file exists", flush=True)
assert path.exists(output_error_checks), "Output checks file does not exist"
assert path.exists(output_error_path), "Output path does not exist"

assert out_error.returncode == 1, "Exit code should be 1"

print(">> Reading json file", flush=True)
with open(output_error_checks, 'r') as f:
    out = json.load(f)
    print(out)


print("All checks succeeded!")
