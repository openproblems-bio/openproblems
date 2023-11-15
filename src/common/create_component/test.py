import os
import subprocess
from os import path
from ruamel.yaml import YAML

## VIASH START
meta = {
    'executable': 'foo'
}
## VIASH END

opv2 = f"{meta['resources_dir']}/openproblems-v2"
output_path = f"{opv2}/src/tasks/label_projection/methods/test_method"

cmd = [
    meta['executable'],
    '--task', 'label_projection',
    '--type', 'method',
    '--name', 'test_method',
    '--language', 'python'
]

print('>> Running the script as test', flush=True)
out = subprocess.run(cmd, stderr=subprocess.STDOUT, cwd=opv2)

if out.stdout:
    print(out.stdout)

if out.returncode:
    print(f"script: '{cmd}' exited with an error.")
    exit(out.returncode)

print('>> Checking whether output files exist', flush=True)
assert os.path.exists(output_path), "Output dir does not exist"

conf_f = path.join(output_path, 'config.vsh.yaml')
assert os.path.exists(conf_f), "Config file does not exist"

script_f = path.join(output_path, "script.py")
assert os.path.exists(script_f), "Script file does not exist"

print('>> Checking file contents', flush=True)
yaml = YAML(typ='safe', pure=True)
with open(conf_f) as f:
    conf_data = yaml.load(f)

assert conf_data['functionality']['name'] == 'test_method', "Name should be equal to 'test_method'"
# assert conf_data['platforms'][0]['image'] == 'python:3.10', "Python image should be equal to python:3.10"


print('All checks succeeded!', flush=True)

