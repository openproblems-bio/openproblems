#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

DATASETS_DIR="resources_test/predict_modality"
OUTPUT_DIR="output/predict_modality"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi
# run benchmark
export NXF_VER=23.04.2

nextflow run . \
  -main-script target/nextflow/predict_modality/workflows/run_benchmark/main.nf \
  -profile docker \
  -resume \
  -entry auto \
  --input_states "$DATASETS_DIR/**/state.yaml" \
  --rename_keys 'input_train_mod1:output_train_mod1,input_train_mod2:output_train_mod2,input_test_mod1:output_test_mod1,input_test_mod2:output_test_mod2' \
  --settings '{"output": "scores.tsv"}' \
  --publish_dir "$OUTPUT_DIR"