#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

OUT_DIR="resources_test/openproblems/task_results_v4"

echo ">>> Fetching raw results..."
aws s3 sync --profile op \
  s3://openproblems-data/resources/task_batch_integration/results/run_2025-01-23_18-03-16/ \
  "$OUT_DIR/raw/" \
  --delete

echo
echo ">>> Processing results..."
if [ -d "$OUT_DIR/processed" ]; then rm -Rf $OUT_DIR/processed; fi
nextflow run target/nextflow/reporting/process_task_results/main.nf \
  -profile docker \
  --input_task_info $OUT_DIR/raw/task_info.yaml \
  --input_dataset_info $OUT_DIR/raw/dataset_uns.yaml \
  --input_method_configs $OUT_DIR/raw/method_configs.yaml \
  --input_metric_configs $OUT_DIR/raw/metric_configs.yaml \
  --input_scores $OUT_DIR/raw/score_uns.yaml \
  --input_trace $OUT_DIR/raw/trace.txt \
  --output_state state.yaml \
  --publishDir $OUT_DIR/processed

echo ">>> Uploading processed results to S3..."
aws s3 sync --profile op \
  "resources_test/openproblems/task_results_v4/" \
  "s3://openproblems-data/resources_test/openproblems/task_results_v4/" \
  --delete --dryrun

echo
echo ">>> Done!"
