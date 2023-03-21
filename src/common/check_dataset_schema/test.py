import subprocess
from os import path
import json

input_path = meta["resources_dir"] + "resources_test/common/pancreas/dataset.h5ad"
input_correct_schema = meta["resources_dir"] +  "resources_test/common/check_schema/anndata_correct.yaml"
input_error_schema = meta["resources_dir"] + "resources_test/common/check_schema/anndata_error.yaml"
output_checks = "checks.json"
output_path = "output.h5ad"


cmd = [
    meta['executable'],
    "--input", input_path,
    "--schema", input_correct_schema,
    "--checks", output_checks,
    "--output", output_path,
]

print(">> Running script as test", flush=True)
subprocess.run(cmd, check=True)

print(">> Checking whether output file exists", flush=True)
assert path.exists(output_checks)
assert path.exists(output_path)

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
    "--checks", output_checks,
    "--output", output_path,
]

print(">> Running script as test", flush=True)
out_error = subprocess.run(cmd_error)

print(">> Checking whether output file exists", flush=True)
assert path.exists(output_checks)
assert path.exists(output_path)

assert out_error.returncode == 1

print(">> Reading json file", flush=True)
with open(output_checks, 'r') as f:
    out = json.load(f)
    print(out)


print("All checks succeeded!")
