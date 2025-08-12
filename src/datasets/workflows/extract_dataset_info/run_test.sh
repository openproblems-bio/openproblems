#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# export TOWER_WORKSPACE_ID=53907369739130

OUTPUT_DIR="output/temp"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

DATASETS_DIR="resources_test/common"

export NXF_VER=22.04.5
nextflow run . \
  -main-script target/nextflow/datasets/workflows/extract_dataset_info/main.nf \
  -profile docker \
  -resume \
  -c src/wf_utils/labels_ci.config \
  -entry auto \
  --input_states "$DATASETS_DIR/**/state.yaml" \
  --rename_keys 'input:output_dataset' \
  --settings '{"output": "dataset_info.yaml"}' \
  --publish_dir "$OUTPUT_DIR" \
  --output_state "state.yaml"