from openproblems.utils import strip_margin
from os import path
import subprocess
import yaml

test_data = strip_margin(f'''\
    |functionality:
    |  name: "phate"
    |  info:
    |    label: PHATE
    |    summary: Preservating trajectories in a dataset by using heat diffusion potential.
    |    description: |
    |      PHATE uses the potential of heat diffusion to preserve trajectories in a dataset via a diffusion process
    |    reference: "moon2019visualizing"
    |    repository_url: "https://github.com/KrishnaswamyLab/PHATE"
    |    documentation_url: "https://github.com/KrishnaswamyLab/PHATE#readme"
    |    preferred_normalization: sqrt_cp10k
    |  # component specific arguments
    |  arguments:
    |    - name: '--n_pca_dims'
    |      type: integer
    |      description: Number of principal components of PCA to use.
    |  resources:
    |    - type: python_script
    |      path: script.py
    |platforms:
    |  - type: docker
    |    image: ghcr.io/openproblems-bio/base_python:1.0.4
    |  - type: nextflow
    |    directives: 
    |      label: [midtime, highmem, highcpu]
    |'''
)

input = "input.vsh.yaml"
with open(input, "w") as file:
    file.write(test_data)

output = "output.vsh.yaml"

cmd = [
    meta['executable'],
    '--input', input,
    '--output', output
]

print('>> Running the script as test', flush=True)
out = subprocess.run(cmd, stderr=subprocess.STDOUT)

if out.returncode:
    print(f"script: '{cmd}' exited with an error.")
    exit(out.returncode)

print('>> Checking whether output files exist', flush=True)
assert path.exists(output), "Output file does not exist"

print('>> Checking file contents', flush=True)
with open(output) as f:
    conf_data = yaml.safe_load(f)

assert "functionality" not in conf_data, ".functionality not removed"
assert "engines" in conf_data, ".platforms not updated"

print('All checks succeeded!', flush=True)
