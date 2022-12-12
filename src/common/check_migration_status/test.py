import subprocess
from os import path
from yaml import load, CSafeLoader

input_sha = meta["resources_dir"] + "temp/openproblems-v1.json"
input_method_info = meta["resources_dir"] +  "temp/method_info.json"
output_path = "output.csv"

cmd = [
    meta['executable'],
    "--git_sha", input_sha,
    "--comp_info", input_method_info,
    "--output", output_path,
]

print(">> Running script as test")
out = subprocess.run(cmd, check=True, capture_output=True, text=True)

print(">> Checking whether output file exists")
assert path.exists(output_path)

print(">> Reading yaml file")
with open(output_path, 'r') as f:
    out = load(f, Loader=CSafeLoader)
    print(out)

print("All checks succeeded!")
