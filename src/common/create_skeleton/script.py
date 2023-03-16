from typing import Any
from ruamel.yaml import YAML
from pathlib import Path
import os
import re

## VIASH START
par = {
  'task': 'denoising',
  'type': 'metric',
  'language': 'python',
  'name': 'new_comp',
  'output': 'src/denoising/methods/new_comp',
  'api_file': 'src/denoising/api/comp_method.yaml'
}
## VIASH END

# TODO: fix adata -> detect input file

def strip_margin(text: str) -> str:
  return re.sub('(\n?)[ \t]*\|', '\\1', text)

def create_config_template(par):
  yaml = YAML()
  config_template = yaml.load(strip_margin(f'''\
    |# The API specifies which type of component this is.
    |# It contains specifications for:
    |#   - The input/output files
    |#   - Common parameters
    |#   - A unit test
    |__merge__: {os.path.relpath(par["api_file"], par["output"])}
    |
    |functionality:
    |  name: {par["name"]}
    |  namespace: {par["task"]}/{par['type']}s
    |
    |  # Metadata for your component (required)
    |  info: 
    |    xx: xx
    |
    |  # Component-specific parameters (optional)
    |  # arguments:
    |  #   - name: "--n_neighbors"
    |  #     type: "integer"
    |  #     default: 5
    |  #     description: Number of neighbors to use.
    |
    |  # Resources required to run the component
    |  resources:
    |    # The script of your component
    |    - type: xx
    |      path: xx
    |    # Additional resources your script needs (optional)
    |    # - type: file
    |    #   path: weights.pt
    |
    |# Target platforms. For more information, see <link to docs>.
    |platforms:
    |  - type: docker
    |    image:
    |    # Add custom dependencies here
    |    setup:
    |  - type: nextflow
    |    directives:
    |      label: [midmem, midcpu]
    |'''
  ))
  return config_template

def add_method_info(conf, par, pretty_name) -> None:
  """Set up the functionality info for a method."""
  conf['functionality']['info'] = {
    'pretty_name': pretty_name,
    'summary': 'FILL IN: A one sentence summary of this method.',
    'description': 'FILL IN: A (multiline) description of how this method works.',
    'reference': 'bibtex_reference_key',
    'documentation_url': 'https://url.to/the/documentation',
    'repository_url': 'https://github.com/organisation/repository',
    'preferred_normalization': 'log_cpm'
  }
  if par["type"] == "control_method":
    del conf['functionality']['info']['reference']
    del conf['functionality']['info']['documentation_url']
    del conf['functionality']['info']['repository_url']

def add_metric_info(conf, par, pretty_name) -> None:
  """Set up the functionality info for a metric."""
  conf['functionality']['info'] = {
    'metrics': [{
      'name': f'{par["name"]}',
      'pretty_name': pretty_name,
      'summary': 'FILL IN: A one sentence summary of this metric.',
      'description': 'FILL IN: A (multiline) description of how this metric works.',
      'reference': 'bibtex_reference_key',
      'documentation_url': 'https://url.to/the/documentation',
      'repository_url': 'https://github.com/organisation/repository',
      'min': 0,
      'max': 1,
      'maximize': 'true',
    }]
  }

def add_script_resource(conf, par) -> None:
  """Add the script to the functionality resources."""
  if par['language'] == 'python':
    conf['functionality']['resources'] = [{
      "type": "python_script",
      "path": "script.py",
    }]
  if par['language'] == 'r':
    conf['functionality']['resources'] = [{
      "type": "r_script",
      "path": "script.R",
    }]

def add_python_setup(conf) -> None:
  """Set up the docker platform for Python."""
  conf['platforms'][0]["image"] = 'python:3.10'
  conf['platforms'][0]["setup"] = [
    {
      "type": "python",
      "pypi": "anndata>=0.8"
    }
  ]

def add_r_setup(conf) -> None:
  """Set up the docker platform for R."""
  conf['platforms'][0]["image"] = 'eddelbuettel/r2u:22.04'
  conf['platforms'][0]["setup"] = [
    {
      "type": "apt",
      "packages": ['libhdf5-dev', 'libgeos-dev', 'python3', 'python3-pip', 'python3-dev', 'python-is-python3']
    },
    {
      "type": "r",
      "cran": "anndata",
      "script": "anndata::install_anndata()" # TODO: does this work?
    }
  ]

def read_api_spec(par):
  with open(par['api_file'], 'r') as f:
    api_spec = YAML().load(f)
  return api_spec

def set_par_values(api_spec) -> dict[str, Any]:
  """Reads in the API file and returns default values for the par object."""
  args = api_spec['functionality']['arguments']
  for argi, arg in enumerate(args):
    arg_type = arg.get("type", "file")
    direction = arg.get("direction", "input")
    key = re.sub("^-*", "", arg['name'])

    # find value
    if arg_type != "file":
      value = arg.get("default", arg.get("example", "..."))
    elif direction == "input":
      key_strip = key.replace("input_", "")
      value = f'resources_test/{par["task"]}/pancreas/{key_strip}.h5ad'
    else:
      key_strip = key.replace("output_", "")
      value = f'{key_strip}.h5ad'

    # store key and value
    api_spec['functionality']['arguments'][argi]["type"] = arg_type
    api_spec['functionality']['arguments'][argi]["direction"] = direction
    api_spec['functionality']['arguments'][argi]["key"] = key
    api_spec['functionality']['arguments'][argi]["value"] = value

def create_python_script(par, api_spec, type):
  args = api_spec['functionality']['arguments']
  par_string = ",\n  ".join(f"'{arg['key']}': '{arg['value']}'" for arg in args)
  read_h5ad_string = "\n".join(
    f"{arg['key']} = ad.read_h5ad(par['{arg['key']}'])"
    for arg in args
    if arg['type'] == "file"
    and arg['direction'] == "input"
  )

  if type == 'metric':
    output_string = strip_margin(f'''\
      |print('Compute metrics', flush=True)
      |# metric_ids and metric_values can have length > 1
      |# but should be of equal length
      |metric_ids = [ '{par['name']}' ]
      |metric_values = [ 0.5 ]
      |
      |print('Create output anndata', flush=True)
      |out = ad.AnnData(
      |    uns={{
      |        'dataset_id': adata.uns['dataset_id'],
      |        'method_id': adata.uns['method_id'],
      |        'metric_ids': metric_ids,
      |        'metric_values': metric_values,
      |    }},
      |)''')
  else:
    # TODO: fix output slots
    output_string = strip_margin('''\
      |print('Preprocess data', flush=True)
      |# ... preprocessing ...
      |
      |print('Train model', flush=True)
      |# ... train model ...
      |
      |print('Generate predictions', flush=True)
      |# ... generate predictions ...
      |
      |print('Create output anndata', flush=True)
      |out = ad.AnnData(
      |    X=y_pred,
      |    uns={
      |        'dataset_id': adata.uns['dataset_id'],
      |        'method_id': meta['functionality_name'],
      |    },
      |)''')

  script = strip_margin(f'''\
    |import anndata as ad
    |
    |## VIASH START
    |par = {{
    |  {par_string}
    |}}
    |meta = {{
    |  'functionality_name': '{par["name"]}'
    |}}
    |## VIASH END
    |
    |print('Reading input files', flush=True)
    |{read_h5ad_string}
    |
    |{output_string}
    |
    |print('Write output to file', flush=True)
    |out.write_h5ad(par['output'], compress='gzip')
    |''')

  return script


def create_r_script(par, api_spec, type):
  args = api_spec['functionality']['arguments']
  par_string = ",\n  ".join(f'{arg["key"]} = "{arg["value"]}"' for arg in args)
  read_h5ad_string = "\n".join(
    f'{arg["key"]} <- anndata::read_h5ad(par[["{arg["key"]}"]])'
    for arg in args
    if arg['type'] == "file"
    and arg['direction'] == "input"
  )

  if type == 'metric':
    output_string = strip_margin(f'''\
      |cat("Compute metrics\\n")
      |# metric_ids and metric_values can have length > 1
      |# but should be of equal length
      |metric_ids <- c( "{par['name']}" )
      |metric_values <- c( 0.5 )
      |
      |cat("Create output anndata\\n")
      |out <- anndata::AnnData(
      |  uns = list(
      |    dataset_id = adata$uns[["dataset_id"]],
      |    method_id = adata$uns[["method_id"]],
      |    metric_ids = metric_ids,
      |    metric_values = metric_values,
      |  )
      |)''')
  else:
    # TODO: fix output slots
    output_string = strip_margin('''\
      |cat("Preprocess data\\n")
      |# ... preprocessing ...
      |
      |cat("Train model\\n")
      |# ... train model ...
      |
      |cat("Generate predictions\\n")
      |# ... generate predictions ...
      |
      |cat("Create output anndata\\n")
      |out <- anndata::AnnData(
      |  X = y_pred,
      |  uns = list(
      |    dataset_id = adata$uns[["dataset_id"]],
      |    method_id = meta[["functionality_name"]],
      |  )
      |)''')

  script = strip_margin(f'''\
    |library(anndata)
    |
    |## VIASH START
    |par <- list(
    |  {par_string}
    |)
    |meta <- list(
    |  functionality_name = "{par["name"]}"
    |)
    |## VIASH END
    |
    |cat("Reading input files\\n")
    |{read_h5ad_string}
    |
    |{output_string}
    |
    |cat("Write output to file\\n")
    |out$write_h5ad(par[["output"]], compress = "gzip")
    |''')

  return script


def main(par):
  ####### CHECK INPUTS #######
  assert re.match("[a-z][a-z0-9_]*", par["name"]), "Name should match the regular expression '[a-z][a-z0-9_]*'. Example: 'my_component'."
  assert len(par['name']) <= 50, "Method name should be at most 50 characters."

  pretty_name = re.sub("_", " ", par['name']).title()

  api_spec = read_api_spec(par)

  ####### CREATE CONFIG #######
  config_template = create_config_template(par)

  # Add component specific info
  if par['type'] == 'metric':
    add_metric_info(config_template, par, pretty_name)
  else:
    add_method_info(config_template, par, pretty_name)

  # add script to resources
  add_script_resource(config_template, par)

  # add elements depending on language
  if par['language'] == 'python':
    add_python_setup(config_template)

  if par['language'] == 'r':
    add_r_setup(config_template)

  ####### CREATE SCRIPT #######
  set_par_values(api_spec)

  if par['language'] == 'python':
    script_out = create_python_script(par, api_spec, par['type'])

  if par['language'] == 'r':
    script_out = create_r_script(par, api_spec, par['type'])

  ####### WRITE OUTPUTS #######
  out_dir = Path(par["output"])
  out_dir.mkdir(exist_ok=True)

  with open(f'{out_dir}/config.vsh.yaml', 'w') as f:
    YAML().dump(config_template, f)

  script_name = config_template['functionality']['resources'][0]['path']

  with open(f'{out_dir}/{script_name}', 'w') as f:
    f.write(script_out)

if __name__ == "__main__":
  main(par)