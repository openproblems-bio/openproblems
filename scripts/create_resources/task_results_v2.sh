#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

OUT_DIR="resources_test/openproblems/task_results_v2"

echo ">> Fetch results in v2 format"
aws s3 sync \
  s3://openproblems-data/resources/batch_integration/results/run_2024-06-28_13-20-27/ \
  "$OUT_DIR/raw/" \
  --delete

echo ">> Process results"
nextflow run openproblems-bio/openproblems \
  -r v2.0.0 \
  -main-script "target/nextflow/common/process_task_results/run/main.nf" \
  -profile docker \
  --input_scores "$OUT_DIR/raw/score_uns.yaml" \
  --input_method_configs "$OUT_DIR/raw/method_configs.yaml" \
  --input_metric_configs "$OUT_DIR/raw/metric_configs.yaml" \
  --input_dataset_info "$OUT_DIR/raw/dataset_uns.yaml" \
  --input_execution "$OUT_DIR/raw/trace.txt" \
  --input_task_info "$OUT_DIR/raw/task_info.yaml" \
  --output_state "state.yaml" \
  --publish_dir "$OUT_DIR/processed/"

echo ">> Uploading results to S3"
aws s3 sync --profile op \
  "resources_test/openproblems/task_results_v2/" \
  "s3://openproblems-data/resources_test/openproblems/task_results_v2/" \
  --delete --dryrun
