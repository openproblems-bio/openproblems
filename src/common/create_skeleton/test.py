import subprocess
from os import path


## VIASH START

meta = {
    'executable': 'foo'
}

## VIASH END

src_path = meta["resources_dir"] + "/openproblems-v2/src"
output_path= 'test_method'


cmd = [
    meta['executable'],
    '--src', src_path,
    '--task', 'label_projection',
    '--comp_type', 'method',
    '--name', 'test_method',
    '--language', 'python',
    '--output', output_path,
]

print('>> Running the script as test', flush=True)
out= subprocess.run(cmd, check=True)

print('>> Checking whether output files exist', flush=True)
assert path.exists(output_path)
conf_f = path.join(output_path, 'config.vsh.yaml')
assert path.exists(conf_f)
script_f = path.join(output_path, "script.py")
assert path.exists(script_f)

print('all checks succeeded!', flush=True)

