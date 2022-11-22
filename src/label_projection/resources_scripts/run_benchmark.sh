#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

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

with open(dataset_dir + "/params_split.yaml", "r") as file:
  split_list = yaml.safe_load(file)
datasets = split_list['param_list']


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
  -resume \
  -params-file "$params_file" \
  --publish_dir "$OUTPUT_DIR"