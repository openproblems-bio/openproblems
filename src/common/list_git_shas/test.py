import subprocess
from os import path

input_path = meta["resources_dir"] + "src/common"
output_path = "output.json"

cmd = [
    meta['executable'],
    "--input_train", input_path,
    "--output", output_path
]

print(">> Running script as test")
out = subprocess.run(cmd, capture_output=True, text=True)


print(">> Checking whether output file exists")
assert path.exists(output_path)

print("All checks succeeded!")