import anndata as ad
import yaml
import shutil
import json

## VIASH START
par = {
  'input': 'resources_test/common/pancreas/dataset.h5ad',
  'schema': 'src/tasks/denoising/api/file_common_dataset.yaml',
  'stop_on_error': False,
  'checks': 'output/error.json',
  'output': 'output/output.h5ad',
  'meta': 'output/meta.json',
}
## VIASH END

def check_structure(slot_info, adata_slot):
  missing = []
  for obj in slot_info:
    if obj['name'] not in adata_slot:
      missing.append(obj['name'])
  return missing

print('Load data', flush=True)
adata = ad.read_h5ad(par['input'])

# create data structure
out = {
  "exit_code": 0,
  "error": {}
}

def is_atomic(obj):
  return isinstance(obj, str) or isinstance(obj, int) or isinstance(obj, bool)

def is_list_of_atomics(obj):
  if not isinstance(obj, list):
    return False
  return all(is_atomic(elem) for elem in obj)

def is_dict_of_atomics(obj):
  if not isinstance(obj, dict):
    return False
  return all(is_atomic(elem) for key, elem in obj.items())


if par['meta'] is not None:
  print("Extract metadata from object", flush=True)
  meta = {
    key: val
    for key, val in adata.uns.items()
    if is_atomic(val) or is_list_of_atomics(val) or is_dict_of_atomics(val)
  }
  with open(par["meta"], "w") as f:
    yaml.dump(meta, f, indent=2)

if par['schema'] is not None:
  print("Check AnnData against schema", flush=True)
  with open(par["schema"], "r") as f:
    data_struct = yaml.safe_load(f)

  def_slots = data_struct['info']['slots']

  out["data_schema"] = "ok"

  for slot in def_slots:
    check = check_structure(def_slots[slot], getattr(adata, slot))
    if bool(check):
      out['exit_code'] = 1
      out['data_schema'] = 'not ok'
      out['error'][slot] = check

  if par['checks'] is not None:
    with open(par["checks"], "w") as f:
      json.dump(out, f, indent=2)

if par['output'] is not None:
  shutil.copyfile(par["input"], par["output"])

if par['stop_on_error']:
  exit(out['exit_code'])  
