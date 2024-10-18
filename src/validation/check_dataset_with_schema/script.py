import anndata as ad
import pandas as pd
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

# TODO: need to refactor to reuse the same helper functions as in 'run_and_check_output.py'.

def check_h5ad_struct(struc, struc_fields, adata_slot):
    missing = []
    if struc == "X":
        struc_fields["name"] = "X"
        struc_fields = [struc_fields]
    for obj in struc_fields:
        adata_data = adata_slot.get(obj['name']) if struc != 'X' else adata_slot
        if obj.get('required') and adata_data is None:
            missing.append(obj['name'])
        # todo: check types
    return missing

def check_df_columns(df, columns):
    missing = []
    for col in columns:
        if col not in df.columns:
            missing.append(col)
    return missing

print("Load schema", flush=True)
with open(par["schema"], "r") as f:
    schema = yaml.safe_load(f)

schema_info = schema.get("info")
assert schema_info, "Schema must contain an 'info' field"

schema_info_format = schema_info.get("format")
assert schema_info_format, "Schema must contain a '.info.format' field"

format_type = schema_info_format.get("type")
assert format_type == "h5ad", ".info.format.type must be 'h5ad'"

# create output data structure
out = {
    "exit_code": 0,
    "error": {},
    "data_schema": "ok"
}

print('Load data', flush=True)
if format_type == "h5ad":
    data = ad.read_h5ad(par['input'])
elif format_type == "csv":
    data = pd.read_csv(par['input'])
elif format_type == "tsv":
    data = pd.read_csv(par['input'], sep="\t")
elif format_type == "parquet":
    data = pd.read_parquet(par['input'])
else:
    raise ValueError(f"Unknown .info.format.type '{format_type}'")

out = {
    "exit_code": 0,
    "error": {},
    "data_schema": "ok"
}
print("Check file against schema", flush=True)
if format_type == "h5ad":
    for struc, struc_fields in schema_info_format.items():
        if struc == "type":
            continue
        print("Checking slot", struc, flush=True)
        missing = check_h5ad_struct(struc, struc_fields, getattr(data, struc))
        if missing:
            print(f"Dataset is missing {struc} {missing}", flush=True)
            out['exit_code'] = 1
            out['data_schema'] = 'not ok'
            out['error'][struc] = missing
elif format_type in ["csv", "tsv", "parquet"]:
    columns = schema_info_format.get("columns") or []
    missing = check_df_columns(data, columns)

with open(par["output"], "w") as f:
    json.dump(out, f, indent=2)

if par['stop_on_error']:
    exit(out['exit_code'])  
