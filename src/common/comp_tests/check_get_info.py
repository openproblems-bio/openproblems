import subprocess
from os import path
import json

## VIASH START
## VIASH END

input_path = meta["resources_dir"] + "/openproblems-v2"
task_id = "denoising"
output_path = "output.json"

cmd = [
    meta['executable'],
    "--input", input_path,
    "--task_id", task_id,
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
assert path.exists(output_path), "Output does not exist"

print(">> Reading json file", flush=True)
with open(output_path, 'r') as f:
    out = json.load(f)
    print(out)

print("All checks succeeded!", flush=True)