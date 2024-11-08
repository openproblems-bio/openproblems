import anndata as ad
import yaml
import numpy as np
import pandas as pd
import scipy
import os
import datetime

## VIASH START
par = {
  'input': 'resources_test/common/pancreas/dataset.h5ad',
  'schema': 'src/datasets/api/file_raw.yaml',
  'output': 'output/meta.yaml',
  'uns_length_cutoff': 100
}
## VIASH END

print('Load data', flush=True)
adata = ad.read_h5ad(par['input']).copy()

if par["schema"]:
  print("Load schema", flush=True)
  with open(par["schema"], "r") as f:
    schema = yaml.safe_load(f)
  
  schema_info = schema.get("info") or {}
  assert schema_info, "Schema must contain an 'info' field"

  schema_info_format = schema_info.get("format") or {}
  assert schema_info_format, "Schema must contain a '.info.format' field"

  assert schema_info_format.get("type") == "h5ad", ".info.format.type must be 'h5ad'"
else:
  schema = None

####################################################################################################
## Helper functions for extracting the dataset metadata in uns                                    ##
####################################################################################################
def is_atomic(obj):
  return pd.api.types.is_scalar(obj)

def to_atomic(obj):
  if isinstance(obj, (np.int32,np.int64)):
    return int(obj)
  elif isinstance(obj, np.float32,np.float64):
    return float(obj)
  elif isinstance(obj, np.bool_):
    return bool(obj)
  elif isinstance(obj, np.str_):
    return str(obj)
  return obj

def is_list_of_atomics(obj):
  if not isinstance(obj, (list, pd.core.series.Series, np.ndarray)):
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


####################################################################################################
## Helper functions for extracting metadata about the used data structures                        ##
####################################################################################################
def get_structure_shape(obj) -> list:
  if isinstance(obj, np.ndarray):
    return list(obj.shape)
  elif scipy.sparse.issparse(obj):
    return list(obj.shape)
  elif isinstance(obj, pd.core.frame.DataFrame):
    return list(obj.shape)
  elif isinstance(obj, pd.core.series.Series):
    return list(obj.shape)
  elif isinstance(obj, list):
    return [len(obj)]
  elif isinstance(obj, dict):
    return [len(obj)]
  elif is_atomic(obj):
    return [1]
  return None

def get_structure_type(obj) -> str:
  # return one of: atomic, dataFrame, vector, dict, denseMatrix, sparseMatrix
  if is_atomic(obj):
    return "atomic"
  elif isinstance(obj, (list,pd.core.series.Series)):
    return "vector"
  elif isinstance(obj, dict):
    return "dict"
  elif isinstance(obj, pd.core.frame.DataFrame):
    return "dataframe"
  elif scipy.sparse.issparse(obj):
    return "sparsematrix"
  elif isinstance(obj, np.ndarray):
    return "densematrix"
  return "other: " + str(type(obj))

def get_structure_dtype(obj) -> str:
  if isinstance(obj, np.ndarray):
    return obj.dtype.name
  elif isinstance(obj, pd.core.series.Series):
    return obj.dtype.name
  elif isinstance(obj, pd.core.frame.DataFrame):
    return [dtype.name for dtype in obj.dtypes]
  elif scipy.sparse.issparse(obj):
    return obj.dtype.name
  elif is_atomic(obj):
    return type(obj).__name__
  return None

def get_structure_schema_info(struct, key) -> dict:
  if schema is None:
    return {}
  
  struct_args = schema_info_format.get(struct, {})
  if struct_args is None:
    return {}
  if struct == "X":
    return struct_args
  
  # look for item with the correct name
  struct_results = [x for x in struct_args if x.get("name") == key]

  # return None if no match is found
  if len(struct_results) != 1:
    return {}

  return struct_results[0]

def get_structure(adata, struct):
  adata_struct = getattr(adata, struct)

  # turn `adata_struct` into a dict for `X`
  if (struct == "X"):
    adata_struct = {"X": adata_struct} if adata_struct is not None else {}

  output = []

  for key, value in adata_struct.items():
    out = {
      "name": key,
      "type": get_structure_type(value),
      "shape": get_structure_shape(value),
      "dtype": get_structure_dtype(value),
    }

    # see if the schema has information about this struct
    schema_info = get_structure_schema_info(struct, key)

    copy = {
      "description": "description",
      "summary": "summary",
      "label": "label",
      "schema_type": "type"
    }
    for k, v in copy.items():
      if schema_info.get(v):
        out[k] = schema_info.get(v)

    output.append(out)
  
  return output

####################################################################################################
## Other helper functions                                                                         ##
####################################################################################################

def get_file_size(path: str) -> int:
  """Get the file size in bytes of the file at the given path."""
  return os.path.getsize(path)

def get_file_creation_time(path: str) -> str:
  """Get the creation time of the file at the given path."""
  # Get file creation time
  creation_time = os.path.getctime(path)
  # Convert creation time from seconds since epoch to a readable timestamp
  creation_time = datetime.datetime.fromtimestamp(creation_time)
  # Format the datetime object as 'DD-MM-YYYY'
  creation_time = creation_time.strftime('%d-%m-%Y')
  return str(creation_time)

print("Extract metadata from object", flush=True)
# Extract metadata about the adata object
uns = {}
for key, val in adata.uns.items():
  if is_atomic(val):
    uns[key] = to_atomic(val)
  elif is_list_of_atomics(val) and len(val) <= par["uns_length_cutoff"]:
    uns[key] = to_list_of_atomics(val)
  elif is_dict_of_atomics(val) and len(val) <= par["uns_length_cutoff"]:
    uns[key] = to_dict_of_atomics(val)

uns["file_size"] = get_file_size(par["input"])
uns["date_created"] = get_file_creation_time(par["input"])

# Extract metadata about the data structures
structure = {
  struct: get_structure(adata, struct)
  for struct
  in ["X", "obs", "var", "obsp", "varp", "obsm", "varm", "layers", "uns"]
}

# ¢reate metadata object
meta = {"uns": uns, "structure": structure}

print("Write metadata to file", flush=True)
with open(par["output"], "w") as f:
  yaml.dump(meta, f, indent=2)
