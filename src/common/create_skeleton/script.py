from ruamel.yaml import YAML
from pathlib import Path


## VIASH START
# The following code has been auto-generated by Viash.
par = {
  'src': './src',
  'task': 'denoising',
  'comp_type': 'metric',
  'language': 'python',
  'name': 'new_comp',
}
meta = {
}

## VIASH END


def add_metric_config(tmpl):
  
  tmpl['functionality']['info']['metrics'] = [{
      'metric_id': 'metric_id',
      'metric_name': 'Metric Name',
      'metric_description': 'metric description',
      'min': 0,
      'max': 1,
      'maximize': 'true',
    }]
  
  return tmpl

def add_method_config(tmpl):

  tmpl['functionality']['info'].update({
    'method_name': 'Method name',
    'preferred_normalization': '',
    'variants': {
      par['name']: '',
      'method_variant1': {
        'preferred_normalization': ''
      }
    }
  })

  return tmpl

def add_python_setup(conf):

  conf['functionality']['resources'][0]['type'] = 'python_script'
  conf['functionality']['resources'][0]['path'] = 'script.py'

  conf['functionality']['test_resources'][0]['type'] = 'python_script'
  conf['functionality']['test_resources'][0]['path'] = 'script.py'

  for i, platform in enumerate(conf['platforms']):
    if platform['type'] == 'docker':
      conf['platforms'][i]['image'] = 'python:3.10'

  return conf

def add_r_setup(conf):

  conf['functionality']['resources'][0]['type'] = 'r_script'
  conf['functionality']['resources'][0]['path'] = 'script.R'

  conf['functionality']['test_resources'][0]['type'] = 'r_script'
  conf['functionality']['test_resources'][0]['path'] = 'script.R'

  for i, platform in enumerate(conf['platforms']):
    if platform['type'] == 'docker':
      pltf = conf['platforms'][i]
      pltf['image'] = 'eddelbuettel/r2u:22.04'
      pltf['setup'].append(
        {
          'type': 'r',
          'cran': [ 'anndata'],
          'bioc': ''
        })
      pltf['setup'].append({
          'type': 'apt',
          'packages': ['libhdf5-dev', 'libgeos-dev', 'python3', 'python3-pip', 'python3-dev', 'python-is-python3']
        })

  return conf


def create_python_script(tmpl_par, comp_type):
  newline = "\n"
  script_templ = f'''import anndata as ad
## VIASH START

par = {{
  # Required arguments for the task
  {newline.join(f"'{key}': '{value}'," for key, value in tmpl_par.items())}
  # Optional method-specific arguments
  'n_neighbors': 5,
  }}

meta = {{
  'functionality_name': 'foo'
}}

## VIASH END

## Data reader
print('Reading input files', flush=True)

adata = ad.read_h5ad(par['{list(templ_par.keys())[0]}'])

print('processing Data', flush=True)
# ... preprocessing ... 
# ... train model ...
# ... generate predictions ...

'''

  if comp_type == 'metric':
    script_templ = script_templ + '''# write output to file
out = ad.AnnData(
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'method_id': adata.uns['method_id'],
        'metric_values': [''],
        'metric_ids': [meta['functionality_name']], # if multiple values, add ids explicitly e.g. ['asw', 'asw_batch']
    },
)

print('writing to output files', flush=True)
out.write_h5ad(par['output'], compress='gzip')
    '''
  else :
    script_templ = script_templ + '''# write output to file
out = ad.AnnData(
    X=y_pred,
    uns={
        'dataset_id': adata.uns['dataset_id'],
        'method_id': meta['functionality_name'],
    },
)

print('writing to output files', flush=True)
out.write_h5ad(par['output'], compress='gzip')
    '''


  return script_templ

def create_r_script(tmpl_par, comp_type):
  newline = "\n"
  script_templ = f'''library(anndata, warn.conflicts = FALSE)

## VIASH START

par <- list(
  # Required arguments for the task
  {newline.join(f'{key} = "{value}",' for key, value in tmpl_par.items())}
  # Optional method-specific arguments
  n_neighbors = 5,
)

meta <- list(
  functionality_name = "foo"
)

## VIASH END

## Data reader
cat("Reading input files\\n")
adata <- read_h5ad(par["{list(templ_par.keys())[0]}"])

cat("processing Data\\n")
# ... preprocessing ...
# ... train model ...
# ... generate predictions ...

'''

  if comp_type == 'metric':
    script_templ = script_templ + '''# write output to file
out <- anndata::AnnData(
  shape = c(0, 0),
  uns = list(
    dataset_id = adata$uns$dataset_id,
    method_id = adata$uns$method_id,
    metric_values = list(""),
    metric_ids = list(meta$functionality_name), # if multiple values, add ids explicitly e.g. list('asw', 'asw_batch')
  )
)

out("writing to output files\\n")
zzz <- adata$write_h5ad(par$output, compression = "gzip")
    '''
  
  else:
        script_templ = script_templ + '''# write output to file
out <- anndata::AnnData(
  X = y_pred,
  uns = list(
    dataset_id = adata$uns$dataset_id,
    method_id = meta$functionality_name,
  )
)

out("writing to output files\\n")
zzz <- adata$write_h5ad(par$output, compression = "gzip")
    '''





  return script_templ


## Create config file
if 'control' in par['comp_type']:
  merge = 'control_method'
else:
  merge = par['comp_type']

config_tmpl = f'''
# points to global config e.g. parameters
__merge__: ../../api/comp_{merge}.yaml
functionality:
  # a unique name for your method, same as what is being output by the script.
  # must match the regex [a-z][a-z0-9_]*
  name: {par['name']}
  namespace: {par["task"]}/{merge}s
  # metadata for your method
  description: A description for your method.
  info: 
    type: {par["comp_type"]}

  # component parameters
  arguments:
    # Method-specific parameters. 
    # Change these to expose parameters of your method to Nextflow (optional)
    - name: "--n_neighbors"
      type: "integer"
      default: 5
      description: Number of neighbors to use.

  # files your script needs
  resources:
    # the script itself
    - type: 
      path:
    # additional resources your script needs (optional)
    - type: file
      path: weights.pt

  # resources for unit testing your component
  test_resources:
    - type: python_script
      path: test.py
    - path: sample_data

# target platforms
platforms:
  # By specifying 'docker' platform, viash will build a standalone
  # executable which uses docker in the back end to run your method.
  - type: docker
    #  you need to specify a base image that contains at least bash and python
    image:
    # You can specify additional dependencies with 'setup'.
    setup:
      - type: python
        pip:
          - pyyaml
          - anndata>=0.8

  # By specifying a 'nextflow', viash will also build a viash module
  # which uses the docker container built above to also be able to
  # run your method as part of a nextflow pipeline.
  - type: nextflow
    directives:
      label: ['midmem', 'midcpu']
'''

yaml = YAML()
conf_tmpl_dict = yaml.load(config_tmpl)

# Add component specific config data

if par['comp_type'] == 'metric':

  config_out = add_metric_config(conf_tmpl_dict)

else:
  
  config_out = add_method_config(conf_tmpl_dict)

  if par['comp_type'] == 'method':
    config_out['functionality']['info']['paper_reference']= ''


# add elements depending on language
if par['language'] == 'python':

  config_out = add_python_setup(config_out)

if par['language'] == 'r':

  config_out = add_r_setup(config_out)


## Create script template

resource_dir = par['src']

task_api = f'{resource_dir}/{par["task"]}/api'
api_conf = f'{task_api}/comp_{merge}.yaml'

with open(api_conf, 'r') as f:
  api_data = yaml.load(f)

args = api_data['functionality']['arguments']

templ_par = {}

for arg in args:
  templ_par[arg['name'].replace('--','')] = ''

if par['language'] == 'python':

  script_out = create_python_script(templ_par, par['comp_type'])

if par['language'] == 'r':

  script_out = create_r_script(templ_par, par['comp_type'])



## Write output
out_dir= Path(par["output"])

out_dir.mkdir(exist_ok=True)

with open(f'{out_dir}/config.vsh.yaml', 'w') as f:
  yaml.dump(config_out, f)

script_f = config_out['functionality']['resources'][0]['path']

with open(f'{out_dir}/{script_f}', 'w') as fpy:
  fpy.write(script_out)