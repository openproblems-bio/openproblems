#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

COMMON_DATASETS="resources/datasets/openproblems_v1"
OUTPUT_DIR="resources/label_projection/datasets/openproblems_v1"

export NXF_VER=22.04.5

nextflow run . \
  -main-script target/nextflow/label_projection/workflows/process_datasets/main.nf \
  -profile docker \
  -entry auto \
  -resume \
  --input_states "$COMMON_DATASETS/**/state.yaml" \
  --rename_keys 'input:output_dataset' \
  --settings '{"output_train": "$id/train.h5ad", "output_test": "$id/test.h5ad", "output_solution": "$id/solution.h5ad"}' \
  --publish_dir "$OUTPUT_DIR" \
  --output_state '$id/state.yaml'
# output_state should be moved to settings once workaround is solved