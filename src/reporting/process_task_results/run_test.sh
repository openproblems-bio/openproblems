#!/bin/bash

# fail on error
set -e

# ensure we're in the root of the repo
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

for TASK in "task_perturbation_prediction"; do
# for TASK in "label_projection"; do
  BASE_DIR="s3://openproblems-data/resources/$TASK/results/"
  
  # find subdir in bucket with latest date
  DATE=$(aws s3 ls $BASE_DIR --recursive | awk '{print $4}' | grep 'task_info.yaml' | sort -r | head -n 1 | sed 's#.*/run_\(.*\)/[^/]*$#\1#')
  
  INPUT_DIR="$BASE_DIR/run_$DATE"
  TASK_STRIP_PREFIX=$(echo $TASK | sed 's/task_//')
  OUTPUT_DIR="../website/results/$TASK_STRIP_PREFIX/data"

  # # temp sync
  # aws s3 sync $INPUT_DIR output/temp

  echo "Processing $TASK - $DATE"

  # start the run
  NXF_VER=23.10.0 nextflow run . \
    -main-script target/nextflow/reporting/process_task_results/main.nf \
    -profile docker \
    -resume \
    -c common/nextflow_helpers/labels_ci.config \
    --id "process" \
    --input_scores "$INPUT_DIR/score_uns.yaml" \
    --input_dataset_info "$INPUT_DIR/dataset_uns.yaml" \
    --input_method_configs "$INPUT_DIR/method_configs.yaml" \
    --input_metric_configs "$INPUT_DIR/metric_configs.yaml" \
    --input_execution "$INPUT_DIR/trace.txt" \
    --input_task_info "$INPUT_DIR/task_info.yaml" \
    --output_state "state.yaml" \
    --publish_dir "$OUTPUT_DIR"

  # cause quarto rerender to index page when in preview mode
  touch ../website/results/$TASK/index.qmd
done