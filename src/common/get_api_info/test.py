import subprocess
from os import path
import json

input_path = meta["resources_dir"] + "src/label_projection"
output_path = "output.json"

cmd = [
    meta['executable'],
    "--input", input_path,
    "--output", output_path,
]

print(">> Running script as test")
out = subprocess.run(cmd, check=True, capture_output=True, text=True)

print(">> Checking whether output file exists")
assert path.exists(output_path)

print(">> Reading json file")
with open(output_path, 'r') as f:
    out = json.load(f)
    print(out)

print("All checks succeeded!")