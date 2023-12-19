#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

DATASETS_DIR="resources_test/predict_modality/neurips2021_bmmc_cite"
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
  -c src/wf_utils/labels_ci.config \
  --input_states "$DATASETS_DIR/state.yaml" \
  --rename_keys 'input_train_mod1:output_train_mod1,input_train_mod2:output_train_mod2,input_test_mod1:output_test_mod1,input_test_mod2:output_test_mod2' \
  --settings '{"output_scores": "scores.yaml", "output_dataset_info": "dataset_info.yaml", "output_method_configs": "method_configs.yaml", "output_metric_configs": "metric_configs.yaml", "output_task_info": "task_info.yaml"}' \
  --publish_dir "$OUTPUT_DIR" \
  --output_state '$id/state.yaml'