#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

DATASETS_DIR="resources_test/spatially_variable_genes"
OUTPUT_DIR="output/temp"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

nextflow run . \
  -main-script target/nextflow/spatially_variable_genes/workflows/run_benchmark/main.nf \
  -profile docker \
  -resume \
  -entry auto \
  -c src/wf_utils/labels_ci.config \
  --input_states "$DATASETS_DIR/**/state.yaml" \
  --rename_keys 'input_dataset:output_dataset,input_solution:output_solution' \
  --publish_dir "$OUTPUT_DIR" \
  --output_state "state.yaml" \
  --settings '{"coord_type_moran_i": "generic", "coord_type_sepal": "grid", "max_neighs_sepal": 4}'