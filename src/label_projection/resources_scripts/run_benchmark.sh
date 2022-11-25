#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

export TOWER_WORKSPACE_ID=53907369739130

DATASETS_DIR="resources/label_projection/datasets/openproblems_v1"
OUTPUT_DIR="resources/label_projection/benchmarks/openproblems_v1"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

params_file="$OUTPUT_DIR/params.yaml"

if [ ! -f $params_file ]; then
  python << HERE
import yaml

dataset_dir = "$DATASETS_DIR"
output_dir = "$OUTPUT_DIR"

# read split datasets yaml
with open(dataset_dir + "/params.yaml", "r") as file:
  split_list = yaml.safe_load(file)
datasets = split_list['param_list']

# figure out where train/test/solution files were stored
param_list = []

for dataset in datasets:
  id = dataset["id"]
  input_train = dataset_dir + "/" + id + ".train.h5ad"
  input_test = dataset_dir + "/" + id + ".test.h5ad"
  input_solution = dataset_dir + "/" + id + ".solution.h5ad"

  obj = {
    'id': id, 
    'dataset_id': dataset["dataset_id"],
    'normalization_id': dataset["normalization_id"],
    'input_train': input_train,
    'input_test': input_test,
    'input_solution': input_solution
  }
  param_list.append(obj)

# write as output file
output = {
  "param_list": param_list,
}

with open(output_dir + "/params.yaml", "w") as file:
  yaml.dump(output, file)
HERE
fi

export NXF_VER=22.04.5
bin/nextflow \
  run . \
  -main-script src/label_projection/workflows/run/main.nf \
  -profile docker \
  -params-file "$params_file" \
  --publish_dir "$OUTPUT_DIR" \
  -with-tower

bin/tools/docker/nextflow/process_log/process_log \
  --output "$OUTPUT_DIR/nextflow_log.tsv"

# bin/viash_build -q label_projection -c '.platforms[.type == "nextflow"].directives.tag := "id: $id, args: $args"'
# bin/viash_build -q label_projection -c '.platforms[.type == "nextflow"].directives.tag := "$id"'

# bin/nextflow run . \
#   -main-script target/nextflow/label_projection/control_methods/majority_vote/main.nf \
#   -profile docker \
#   --input_train resources_test/label_projection/pancreas/train.h5ad \
#   --input_test resources_test/label_projection/pancreas/test.h5ad \
#   --input_solution resources_test/label_projection/pancreas/solution.h5ad \
#   --publish_dir foo