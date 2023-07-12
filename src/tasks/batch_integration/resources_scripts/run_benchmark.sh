#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

export TOWER_WORKSPACE_ID=53907369739130

DATASETS_DIR="resources/batch_integration/datasets/openproblems_v1"
OUTPUT_DIR="resources/batch_integration/benchmarks/openproblems_v1"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

params_file="$OUTPUT_DIR/params.yaml"

if [ ! -f $params_file ]; then
  python << HERE
import anndata as ad
import glob
import yaml

h5ad_files = glob.glob("$DATASETS_DIR/**/*.h5ad", recursive=True)

# figure out where dataset files are stored
param_list = []

for h5ad_file in h5ad_files:
  print(f"Checking {h5ad_file}")
  adata = ad.read_h5ad(h5ad_file, backed=True)

  dataset_id = adata.uns["dataset_id"].replace("/", ".")
  normalization_id = adata.uns["normalization_id"]
  id = dataset_id + "." + normalization_id

  obj = {
    'id': id,
    'input': h5ad_file, 
    # 'dataset_id': dataset_id,
    # 'normalization_id': normalization_id
  }
  param_list.append(obj)

# write as output file
output = {
  "param_list": param_list,
}

with open("$params_file", "w") as file:
  yaml.dump(output, file)
HERE
fi

export NXF_VER=22.04.5
nextflow \
  run . \
  -main-script src/tasks/batch_integration/workflows/run/main.nf \
  -profile docker \
  -resume \
  -params-file "$params_file" \
  --publish_dir "$OUTPUT_DIR"
