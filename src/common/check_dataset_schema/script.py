import anndata as ad
import yaml
import shutil
import json
import numpy as np
import pandas as pd

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
    if obj.get('required') and obj['name'] not in adata_slot:
      missing.append(obj['name'])
  return missing

print('Load data', flush=True)
adata = ad.read_h5ad(par['input']).copy()

# create data structure
out = {
  "exit_code": 0,
  "error": {},
  "data_schema": "ok"
}

def is_atomic(obj):
  return isinstance(obj, str) or isinstance(obj, int) or isinstance(obj, bool) or isinstance(obj, float)

def to_atomic(obj):
  if isinstance(obj, np.float64):
    return float(obj)
  elif isinstance(obj, np.int64):
    return int(obj)
  elif isinstance(obj, np.bool_):
    return bool(obj)
  elif isinstance(obj, np.str_):
    return str(obj)
  return obj

def is_list_of_atomics(obj):
  if not isinstance(obj, (list,pd.core.series.Series,np.ndarray)):
    return False
  return all(is_atomic(elem) for elem in obj)

def to_list_of_atomics(obj):
  if isinstance(obj, pd.core.series.Series):
    obj = obj.to_numpy()
  if isinstance(obj, np.ndarray):
    obj = obj.tolist()
  return [to_atomic(elem) for elem in obj]

def is_dict_of_atomics(obj):
  if not isinstance(obj, dict):
    return False
  return all(is_atomic(elem) for _, elem in obj.items())

def to_dict_of_atomics(obj):
  return {k: to_atomic(v) for k, v in obj.items()}

if par['meta'] is not None:
  print("Extract metadata from object", flush=True)
  uns = {}
  for key, val in adata.uns.items():
    if is_atomic(val):
      uns[key] = to_atomic(val)
    elif is_list_of_atomics(val):
      uns[key] = to_list_of_atomics(val)
    elif is_dict_of_atomics(val):
      uns[key] = to_dict_of_atomics(val)
  structure = {
    struct: list(getattr(adata, struct).keys())
    for struct
    in ["obs", "var", "obsp", "varp", "obsm", "varm", "layers", "uns"]
  }
  meta = {"uns": uns, "structure": structure}
  with open(par["meta"], "w") as f:
    yaml.dump(meta, f, indent=2)

if par['schema'] is not None:
  print("Check AnnData against schema", flush=True)
  with open(par["schema"], "r") as f:
    data_struct = yaml.safe_load(f)

  def_slots = data_struct['info']['slots']

  missing= []
  for slot in def_slots:
    missing_x = False
    if slot == "X":
      if adata.X is None:
        missing_x = True
      continue
    missing = check_structure(def_slots[slot], getattr(adata, slot))
    if missing_x:
      missing.append("X")
    if missing:
      out['exit_code'] = 1
      out['data_schema'] = 'not ok'
      out['error'][slot] = missing

  if par['checks'] is not None:
    with open(par["checks"], "w") as f:
      json.dump(out, f, indent=2)

if par['output'] is not None and out["data_schema"] == "ok":
  shutil.copyfile(par["input"], par["output"])

if par['stop_on_error']:
  exit(out['exit_code'])  
