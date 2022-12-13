import subprocess
from os import path
import json

input_sha = meta["resources_dir"] + "resources_test/common/input_git_sha.json"
input_method_info = meta["resources_dir"] +  "resources_test/common/method_info.json"
output_path = "output.json"

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

print(">> Reading json file")
with open(output_path, 'r') as f:
    out = json.load(f)
    print(out)

print("All checks succeeded!")
