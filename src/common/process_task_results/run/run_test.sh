#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# export TOWER_WORKSPACE_ID=53907369739130

DATASETS_DIR="resources/batch_integration/results"
OUTPUT_DIR="../website/results/batch_integration_feature/data"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

export NXF_VER=22.04.5

nextflow run . \
  -main-script target/nextflow/common/workflows/transform_meta/main.nf \
  -profile docker \
  -resume \
  -c src/wf_utils/labels_ci.config \
  --id "get_results_test" \
  --input_scores "$DATASETS_DIR/scores.yaml" \
  --input_dataset_info "$DATASETS_DIR/dataset_info.yaml" \
  --input_method_configs "$DATASETS_DIR/method_configs.yaml" \
  --input_metric_configs "$DATASETS_DIR/metric_configs.yaml" \
  --input_execution "$DATASETS_DIR/trace.txt" \
  --input_task_info "$DATASETS_DIR/task_info.yaml" \
  --task_id "batch_integration" \
  --output_scores "results.json"\
  --output_method_info "method_info.json"\
  --output_metric_info "metric_info.json"\
  --output_dataset_info "dataset_info.json"\
  --output_task_info "task_info.json" \
  --publish_dir "$OUTPUT_DIR"


# nextflow run . \
#   -main-script target/nextflow/common/workflows/transform_meta/main.nf \
#   -profile docker \
#   -resume \
#   -entry auto \
#   --input_states "$DATASETS_DIR/state.yaml" \
#   --rename_keys 'input_scores:output_scores,input_dataset_info:output_dataset_info, input_method_configs:output_method_configs, input_metric_configs:output_metric_configs, ' \
#   --settings '{"task_id": "batch_integration", "output_scores": "results.json", "output_method_info": "method_info.json", "output_metric_info": "metric_info.json", "output_dataset_info": "dataset_info.json", "output_task_info":"task_info.json"}' \
#   --publish_dir "$OUTPUT_DIR"