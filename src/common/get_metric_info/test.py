import subprocess
from os import path
from yaml import load, CSafeLoader

input_path = meta["resources_dir"] + "src/label_projection"
output_path = "output.yaml"

cmd = [
    meta['executable'],
    "--input", input_path,
    "--output", output_path,
]

print(">> Running script as test")
out = subprocess.run(cmd, check=True, capture_output=True, text=True)

print(">> Checking whether output file exists")
assert path.exists(output_path)

print(">> Reading yaml file")
with open(output_path, 'r') as f:
    out = load(f, Loader= CSafeLoader)
    print(out[0])

print("All checks succeeded!")
