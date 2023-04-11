from typing import Any
from ruamel.yaml import YAML
from pathlib import Path
import os
import re
import subprocess

yaml = YAML()
yaml.indent(mapping=2, sequence=4, offset=2)

## VIASH START
par = {
  'task': 'denoising',
  'type': 'method',
  'language': 'python',
  'name': 'new_comp',
  'output': 'src/denoising/methods/new_comp',
  'api_file': 'src/denoising/api/comp_method.yaml'
}
## VIASH END

def strip_margin(text: str) -> str:
  return re.sub('(\n?)[ \t]*\|', '\\1', text)

def create_config_template(par):
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
      "pypi": "anndata~=0.8"
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
      "type": "python",
      "pypi": "anndata~=0.8"
    },
    {
      "type": "r",
      "cran": "anndata"
    }
  ]

def set_par_values(config) -> dict[str, Any]:
  """Adds values to each of the arguments in a config file."""
  args = config['functionality']['arguments']
  for argi, arg in enumerate(args):
    key = re.sub("^-*", "", arg['name'])

    # find value
    if arg["type"] != "file":
      value = arg.get("default", arg.get("example", "..."))
    elif arg["direction"] == "input":
      key_strip = key.replace("input_", "")
      value = f'resources_test/{par["task"]}/pancreas/{key_strip}.h5ad'
    else:
      key_strip = key.replace("output_", "")
      value = f'{key_strip}.h5ad'

    # store key and value
    config['functionality']['arguments'][argi]["key"] = key
    config['functionality']['arguments'][argi]["value"] = value
  
def look_for_adata_arg(args, uns_field):
  """Look for an argument that has a .uns[uns_field] in its info.slots."""
  for arg in args:
    uns = arg.get("info", {}).get("slots", {}).get("uns", [])
    for unval in uns:
      if unval.get("name") == uns_field:
        return arg["key"]
  return "adata"

def write_output_python(arg, copy_from_adata, is_metric):
  """Create code for writing the output h5ad files."""
  slots = arg.get("info", {}).get("slots", {})
  outer = []
  for group_name, slots in slots.items():
    inner = []
    for slot in slots:
      if group_name == "uns" and slot["name"] in ["dataset_id", "normalization_id"]:
        value = f"{copy_from_adata}.uns['{slot['name']}']"
      elif group_name == "uns" and slot["name"] == "method_id":
        if is_metric:
          value = f"{copy_from_adata}.uns['{slot['name']}']"
        else:
          value = "meta['functionality_name']"
      else:
        value = group_name + "_" + slot["name"]
      inner.append(f"'{slot['name']}': {value}")
    inner_values = ',\n    '.join(inner)
    outer.append(f"{group_name}={{\n    {inner_values}\n  }}")
  outer_values = ',\n  '.join(outer)
  return strip_margin(
    f'''\
      |print("Write {arg["key"]} AnnData to file", flush=True)
      |{arg["key"]} = ad.AnnData(
      |  {outer_values}
      |)
      |{arg["key"]}.write_h5ad(par['{arg["key"]}'], compression='gzip')'''
  )

def write_output_r(arg, copy_from_adata, is_metric):
  """Create code for writing the output h5ad files."""
  slots = arg.get("info", {}).get("slots", {})
  outer = []
  for group_name, slots in slots.items():
    inner = []
    for slot in slots:
      if group_name == "uns" and slot["name"] in ["dataset_id", "normalization_id"]:
        value = f"{copy_from_adata}$uns[[\"{slot['name']}\"]]"
      elif group_name == "uns" and slot["name"] == "method_id":
        if is_metric:
          value = f"{copy_from_adata}$uns[[\"{slot['name']}\"]]"
        else:
          value = "meta[[\"functionality_name\"]]"
      else:
        value = group_name + "_" + slot["name"]
      inner.append(f"{slot['name']} = {value}")
    inner_values = ',\n    '.join(inner)
    outer.append(f"{group_name} = list(\n    {inner_values}\n  )")
  outer_values = ',\n  '.join(outer)
  return strip_margin(
    f'''\
      |cat("Write {arg["key"]} AnnData to file\\n")
      |{arg["key"]} <- anndata::AnnData(
      |  {outer_values}
      |)
      |{arg["key"]}$write_h5ad(par[["{arg["key"]}"]], compression = "gzip")'''
  )

def create_python_script(par, config, type):
  args = config['functionality']['arguments']

  # create the arguments of the par string
  par_string = ",\n  ".join(f"'{arg['key']}': '{arg['value']}'" for arg in args)

  # create code for reading the input h5ad file
  read_h5ad_string = "\n".join(
    f"{arg['key']} = ad.read_h5ad(par['{arg['key']}'])"
    for arg in args
    if arg['type'] == "file"
    and arg['direction'] == "input"
  )

  # determine which adata to copy from
  copy_from_adata = look_for_adata_arg(args, "method_id" if type == "metric" else "dataset_id")

  # create code for writing the output h5ad files
  write_h5ad_string = "\n".join(
    write_output_python(arg, copy_from_adata, type == "metric")
    for arg in args
    if arg["type"] == "file"
    and arg["direction"] == "output"
  )

  if type == 'metric':
    processing_string = strip_margin(f'''\
      |print('Compute metrics', flush=True)
      |# metric_ids and metric_values can have length > 1
      |# but should be of equal length
      |uns_metric_ids = [ '{par['name']}' ]
      |uns_metric_values = [ 0.5 ]''')
  else:
    processing_string = strip_margin(f'''\
      |print('Preprocess data', flush=True)
      |# ... preprocessing ...
      |
      |print('Train model', flush=True)
      |# ... train model ...
      |
      |print('Generate predictions', flush=True)
      |# ... generate predictions ...''')

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
    |{processing_string}
    |
    |{write_h5ad_string}
    |''')

  return script

def create_r_script(par, api_spec, type):
  args = api_spec['functionality']['arguments']

  # create the arguments of the par string
  par_string = ",\n  ".join(f'{arg["key"]} = "{arg["value"]}"' for arg in args)

  # create helpers for reading the h5ad file
  read_h5ad_string = "\n".join(
    f'{arg["key"]} <- anndata::read_h5ad(par[["{arg["key"]}"]])'
    for arg in args
    if arg['type'] == "file"
    and arg['direction'] == "input"
  )

  # determine which adata to copy from
  copy_from_adata = look_for_adata_arg(args, "method_id" if type == "metric" else "dataset_id")

  # create code for writing the output h5ad files
  write_h5ad_string = "\n".join(
    write_output_r(arg, copy_from_adata, type == "metric")
    for arg in args
    if arg["type"] == "file"
    and arg["direction"] == "output"
  )

  if type == 'metric':
    processing_string = strip_margin(f'''\
      |cat("Compute metrics\\n")
      |# metric_ids and metric_values can have length > 1
      |# but should be of equal length
      |uns_metric_ids <- c("{par['name']}")
      |uns_metric_values <- c(0.5)''')
  else:
    processing_string = strip_margin(f'''\
      |cat("Preprocess data\\n")
      |# ... preprocessing ...
      |
      |cat("Train model\\n")
      |# ... train model ...
      |
      |cat("Generate predictions\\n")
      |# ... generate predictions ...''')

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
    |{processing_string}
    |
    |{write_h5ad_string}
    |''')

  return script

def read_viash_config(file):
  # read in config
  command = ["viash", "config", "view", str(file)]

  # Execute the command and capture the output
  output = subprocess.check_output(command, universal_newlines=True)

  # Parse the output as YAML
  config = yaml.load(output)

  return config


def main(par):
  ####### CHECK INPUTS #######
  assert re.match("[a-z][a-z0-9_]*", par["name"]), "Name should match the regular expression '[a-z][a-z0-9_]*'. Example: 'my_component'."
  assert len(par['name']) <= 50, "Method name should be at most 50 characters."

  pretty_name = re.sub("_", " ", par['name']).title()

  ####### CREATE OUTPUT DIR #######
  out_dir = Path(par["output"])
  out_dir.mkdir(exist_ok=True)

  ####### CREATE CONFIG #######
  config_file = out_dir / 'config.vsh.yaml'

  # get config template
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

  with open(config_file, 'w') as f:
    yaml.dump(config_template, f)

  ####### CREATE SCRIPT #######
  script_file = out_dir / config_template['functionality']['resources'][0]['path']

  # touch file
  script_file.touch()

  # read config with viash
  final_config = read_viash_config(config_file)

  # set reasonable values
  set_par_values(final_config)

  if par['language'] == 'python':
    script_out = create_python_script(par, final_config, par['type'])

  if par['language'] == 'r':
    script_out = create_r_script(par, final_config, par['type'])
  
  # write script
  with open(script_file, 'w') as f:
    f.write(script_out)


if __name__ == "__main__":
  main(par)