#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# export TOWER_WORKSPACE_ID=53907369739130

DATASETS_DIR="resources/label_projection"
OUTPUT_DIR="output/test"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

export NXF_VER=22.04.5
nextflow run . \
  -main-script target/nextflow/label_projection/workflows/run_benchmark/main.nf \
  -profile docker \
  -resume \
  -entry auto \
  -c src/wf_utils/labels_ci.config \
  --input_states "$DATASETS_DIR/**/state.yaml" \
  --rename_keys 'input_train:output_train,input_test:output_test,input_solution:output_solution' \
  --settings '{"output": "scores.tsv"}' \
  --publish_dir "$OUTPUT_DIR"\
  --output_state '$id/state.yaml'