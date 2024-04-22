from typing import Any
from pathlib import Path
import sys
import os
import re

## VIASH START
par = {
  "task": "denoising",
  "type": "method",
  "language": "python",
  "name": "new_comp",
  "output": "src/tasks/denoising/methods/new_comp",
  "api_file": "src/tasks/denoising/api/comp_method.yaml",
  "viash_yaml": "_viash.yaml"
}
## VIASH END

# import helper function
sys.path.append(meta["resources_dir"])
from read_and_merge_yaml import read_and_merge_yaml

def strip_margin(text: str) -> str:
  return re.sub("(^|\n)[ \t]*\|", "\\1", text)

def create_config(par, component_type, pretty_name, script_path) -> str:
  info_str = generate_info(par, component_type, pretty_name)
  resources_str = generate_resources(par, script_path)
  docker_platform = generate_docker_platform(par)

  return strip_margin(f'''\
    |# The API specifies which type of component this is.
    |# It contains specifications for:
    |#   - The input/output files
    |#   - Common parameters
    |#   - A unit test
    |__merge__: {os.path.relpath(par["api_file"], par["output"])}
    |
    |functionality:
    |  # A unique identifier for your component (required).
    |  # Can contain only lowercase letters or underscores.
    |  name: {par["name"]}
    |
    |  # Metadata for your component
    |  info:
    |{info_str}
    |  # Component-specific parameters (optional)
    |  # arguments:
    |  #   - name: "--n_neighbors"
    |  #     type: "integer"
    |  #     default: 5
    |  #     description: Number of neighbors to use.
    |
    |  # Resources required to run the component
    |  resources:
    |{resources_str}
    |platforms:
    |  # Specifications for the Docker image for this component.
    |{docker_platform}
    |  # This platform allows running the component natively
    |  - type: native
    |  # Allows turning the component into a Nextflow module / pipeline.
    |  - type: nextflow
    |    directives:
    |      label: [midtime,midmem, midcpu]
    |'''
  )

def generate_info(par, component_type, pretty_name) -> str:
  """Generate the functionality info for a component."""
  if component_type in ["method", "control_method"]:
    str = strip_margin(f'''\
      |    # A relatively short label, used when rendering visualisarions (required)
      |    label: {pretty_name}
      |    # A one sentence summary of how this method works (required). Used when 
      |    # rendering summary tables.
      |    summary: "FILL IN: A one sentence summary of this method."
      |    # A multi-line description of how this component works (required). Used
      |    # when rendering reference documentation.
      |    description: |
      |      FILL IN: A (multi-line) description of how this method works.
      |    # Which normalisation method this component prefers to use (required).
      |    preferred_normalization: log_cp10k
      |''')
    if component_type == "method":
      str += strip_margin(f'''\
        |    # A reference key from the bibtex library at src/common/library.bib (required).
        |    reference: bibtex_reference_key
        |    # URL to the documentation for this method (required).
        |    documentation_url: https://url.to/the/documentation
        |    # URL to the code repository for this method (required).
        |    repository_url: https://github.com/organisation/repository
        |''')
    return str
  elif component_type == "metric":
    return strip_margin(f'''\
      |    metrics:
      |      # A unique identifier for your metric (required).
      |      # Can contain only lowercase letters or underscores.
      |      name: {par["name"]}
      |      # A relatively short label, used when rendering visualisarions (required)
      |      label: {pretty_name}
      |      # A one sentence summary of how this metric works (required). Used when 
      |      # rendering summary tables.
      |      summary: "FILL IN: A one sentence summary of this metric."
      |      # A multi-line description of how this component works (required). Used
      |      # when rendering reference documentation.
      |      description: |
      |        FILL IN: A (multi-line) description of how this metric works.
      |      # A reference key from the bibtex library at src/common/library.bib (required).
      |      reference: bibtex_reference_key
      |      # URL to the documentation for this metric (required).
      |      documentation_url: https://url.to/the/documentation
      |      # URL to the code repository for this metric (required).
      |      repository_url: https://github.com/organisation/repository
      |      # The minimum possible value for this metric (required)
      |      min: 0
      |      # The maximum possible value for this metric (required)
      |      max: 1
      |      # Whether a higher value represents a 'better' solution (required)
      |      maximize: true
      |''')


def generate_resources(par, script_path) -> str:
  """Add the script to the functionality resources."""
  if par["language"] == "python":
    type_str = "python_script"
  elif par["language"] == "r":
    type_str = "r_script"

  return strip_margin(f'''\
    |    # The script of your component (required)
    |    - type: {type_str}
    |      path: {script_path}
    |    # Additional resources your script needs (optional)
    |    # - type: file
    |    #   path: weights.pt
    |''')

def generate_docker_platform(par) -> str:
  """Set up the docker platform for Python."""
  if par["language"] == "python":
    image_str = "ghcr.io/openproblems-bio/base_python:1.0.4"
    setup_type = "python"
    package_example = "scib==1.1.5"
  elif par["language"] == "r":
    image_str = "ghcr.io/openproblems-bio/base_r:1.0.4"
    setup_type = "r"
    package_example = "tidyverse"
  return strip_margin(f'''\
    |  - type: docker
    |    image: {image_str}
    |    # Add custom dependencies here (optional). For more information, see
    |    # https://viash.io/reference/config/platforms/docker/#setup .
    |    # setup:
    |    #   - type: {setup_type}
    |    #     packages: {package_example}
    |''')

def set_par_values(config) -> None:
  """Adds values to each of the arguments in a config file."""
  args = config['functionality']['arguments']
  for argi, arg in enumerate(args):
    key = re.sub("^-*", "", arg['name'])

    # find value
    if arg["type"] != "file":
      value = arg.get("default", arg.get("example", "..."))
    elif arg.get("direction", "input") == "input":
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
    and arg.get('direction', "input") == "input"
  )

  # determine which adata to copy from
  copy_from_adata = look_for_adata_arg(args, "method_id" if type == "metric" else "dataset_id")

  # create code for writing the output h5ad files
  write_h5ad_string = "\n".join(
    write_output_python(arg, copy_from_adata, type == "metric")
    for arg in args
    if arg["type"] == "file"
    and arg.get("direction", "input") == "output"
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
    |# Note: this section is auto-generated by viash at runtime. To edit it, make changes
    |# in config.vsh.yaml and then run `viash config inject config.vsh.yaml`.
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
    and arg.get("direction", "input") == "input"
  )

  # determine which adata to copy from
  copy_from_adata = look_for_adata_arg(args, "method_id" if type == "metric" else "dataset_id")

  # create code for writing the output h5ad files
  write_h5ad_string = "\n".join(
    write_output_r(arg, copy_from_adata, type == "metric")
    for arg in args
    if arg["type"] == "file"
    and arg.get("direction", "input") == "output"
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

# def read_viash_config(file):
#   file = file.absolute()

#   # read in config
#   command = ["viash", "config", "view", str(file)]

#   # Execute the command and capture the output
#   output = subprocess.check_output(
#     command,
#     universal_newlines=True,
#     cwd=str(file.parent)
#   )

#   # Parse the output as YAML
#   config = yaml.load(output)

#   return config


def main(par):
  ####### CHECK INPUTS #######
  print("Check inputs", flush=True)
  assert re.match("[a-z][a-z0-9_]*", par["name"]), "Name should match the regular expression '[a-z][a-z0-9_]*'. Example: 'my_component'."
  assert len(par['name']) <= 50, "Method name should be at most 50 characters."

  pretty_name = re.sub("_", " ", par['name']).title()

  ####### CHECK LANGUAGE #######
  print("Check language", flush=True)
  # check language and determine script path
  if par["language"] == "python":
    script_path = "script.py"
  elif par["language"] == "r":
    script_path = "script.R"
  else:
    sys.exit(f"Unrecognized language parameter '{par['language']}'.")

  ## CHECK API FILE
  print("Check API file", flush=True)
  api_file = Path(par["api_file"])
  viash_yaml = Path(par["viash_yaml"])
  project_dir = viash_yaml.parent
  if not api_file.exists():
    comp_types = [x.with_suffix("").name.removeprefix("comp_") for x in api_file.parent.glob("**/comp_*.y*ml")]
    list.sort(comp_types)
    sys.exit(strip_margin(f"""\
      |Error: Invalid --type argument.
      |  Reason: Could not find API file at '{api_file.relative_to(project_dir)}'.
      |  Possible values for --type: {', '.join(comp_types)}."""))
  
  ## READ API FILE
  print("Read API file", flush=True)
  api = read_and_merge_yaml(api_file)
  comp_type = api.get("functionality", {}).get("info", {}).get("type", {})
  if not comp_type:
    sys.exit(strip_margin(f"""\
      |Error: API file is incorrectly formatted.
      |  Reason: Could not find component type at `.functionality.info.type`.'
      |  Please fix the formatting of the API file."""))

  ####### CREATE OUTPUT DIR #######
  print("Create output dir", flush=True)
  out_dir = Path(par["output"])
  out_dir.mkdir(exist_ok=True)

  ####### CREATE CONFIG #######
  print("Create config", flush=True)
  config_file = out_dir / "config.vsh.yaml"

  # get config template
  config_str = create_config(par, comp_type, pretty_name, script_path)

  with open(config_file, "w") as f:
    f.write(config_str)

  ####### CREATE SCRIPT #######
  print("Create script", flush=True)
  script_file = out_dir / script_path

  # set reasonable values
  set_par_values(api)

  if par["language"] == "python":
    script_out = create_python_script(par, api, comp_type)

  if par["language"] == "r":
    script_out = create_r_script(par, api, comp_type)
  
  # write script
  with open(script_file, "w") as f:
    f.write(script_out)

  print("Done!", flush=True)


if __name__ == "__main__":
  main(par)
