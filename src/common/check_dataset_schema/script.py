import anndata as ad
import yaml
import json

## VIASH START
par = {
  'input': 'work/d4/f4fabc8aa4f2308841d4ab57bcff62/_viash_par/input_1/dataset.h5ad',
  'schema': 'work/d4/f4fabc8aa4f2308841d4ab57bcff62/_viash_par/schema_1/schema.yaml',
  'stop_on_error': False,
  'output': 'work/d4/f4fabc8aa4f2308841d4ab57bcff62/out.yaml',
}
## VIASH END

def check_structure(slot, slot_info, adata_slot):
  missing = []
  if slot == "X":
    slot_info["name"] = "X"
    slot_info = [slot_info]
  for obj in slot_info:
    adata_data = adata_slot.get(obj['name']) if slot != 'X' else adata_slot
    if obj.get('required') and adata_data is None:
      missing.append(obj['name'])
    # todo: check types
  return missing

print('Load data', flush=True)
adata = ad.read_h5ad(par['input'])

# create data structure
out = {
  "exit_code": 0,
  "error": {},
  "data_schema": "ok"
}

print("Check AnnData against schema", flush=True)
with open(par["schema"], "r") as f:
  data_struct = yaml.safe_load(f)

def_slots = data_struct['info']['slots']

out = {
  "exit_code": 0,
  "error": {},
  "data_schema": "ok"
}
for slot in def_slots:
  print("Checking slot", slot, flush=True)
  missing = check_structure(slot, def_slots[slot], getattr(adata, slot))
  if missing:
    out['exit_code'] = 1
    out['data_schema'] = 'not ok'
    out['error'][slot] = missing

with open(par["output"], "w") as f:
  json.dump(out, f, indent=2)

if par['stop_on_error']:
  exit(out['exit_code'])  
