#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

OUT_DIR="resources_test/openproblems/task_results_v4"
RUN_S3_PATH="s3://openproblems-data/resources/task_batch_integration/results/run_2025-01-23_18-03-16"

# Derive an ISO 8601 timestamp from the run directory name
# (run_YYYY-MM-DD_HH-MM-SS -> YYYY-MM-DDTHH:MM:SSZ). The raw task_info.yaml does
# not carry a timestamp, so it must be passed explicitly to process_task_results.
RUN_DATE=$(basename "$RUN_S3_PATH")   # run_2025-01-23_18-03-16
RUN_DATE=${RUN_DATE#run_}             # 2025-01-23_18-03-16
DATE_PART=${RUN_DATE%%_*}             # 2025-01-23
TIME_PART=${RUN_DATE#*_}              # 18-03-16
TIME_PART=${TIME_PART//-/:}           # 18:03:16
TIMESTAMP="${DATE_PART}T${TIME_PART}Z"

echo ">>> Fetching raw results..."
aws s3 sync --profile op \
  "$RUN_S3_PATH/" \
  "$OUT_DIR/raw/" \
  --delete

echo
echo ">>> Processing results (timestamp: $TIMESTAMP)..."
if [ -d "$OUT_DIR/processed" ]; then rm -Rf $OUT_DIR/processed; fi
nextflow run target/nextflow/reporting/process_task_results/main.nf \
  -profile docker \
  --input_task_info $OUT_DIR/raw/task_info.yaml \
  --input_dataset_info $OUT_DIR/raw/dataset_uns.yaml \
  --input_method_configs $OUT_DIR/raw/method_configs.yaml \
  --input_metric_configs $OUT_DIR/raw/metric_configs.yaml \
  --input_scores $OUT_DIR/raw/score_uns.yaml \
  --input_trace $OUT_DIR/raw/trace.txt \
  --timestamp "$TIMESTAMP" \
  --output_state state.yaml \
  --publishDir $OUT_DIR/processed

echo ">>> Uploading processed results to S3..."
aws s3 sync --profile op \
  "resources_test/openproblems/task_results_v4/" \
  "s3://openproblems-data/resources_test/openproblems/task_results_v4/" \
  --delete --dryrun

echo
echo ">>> Done!"
