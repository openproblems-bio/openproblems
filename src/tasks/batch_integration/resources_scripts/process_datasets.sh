#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

COMMON_DATASETS="resources/datasets/openproblems_v1"
OUTPUT_DIR="resources/batch_integration/datasets/openproblems_v1"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

export NXF_VER=22.04.5
nextflow run . \
  -main-script src/tasks/batch_integration/workflows/process_datasets/main.nf \
  -profile docker \
  -entry auto \
  -resume \
  --id resources \
  --input_states "resources/datasets/openproblems_v1/**/state.yaml" \
  --rename_keys 'input:output_dataset' \
  --settings '{"output_dataset": "dataset.h5ad", "output_solution": "solution.h5ad"}' \
  --publish_dir "$OUTPUT_DIR"