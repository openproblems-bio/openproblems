#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

DATASETS_DIR="resources_test/batch_integration"
OUTPUT_DIR="resources_test/common/task_metadata"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

echo ">> Running benchmark"
# TODO: manually generate contents of output_dir
export NXF_VER=22.04.5

nextflow run . \
  -main-script target/nextflow/batch_integration/workflows/run_benchmark/main.nf \
  -profile docker \
  -resume \
  -c src/wf_utils/labels_ci.config \
  -with-trace \
  -entry auto \
  --input_states "$DATASETS_DIR/pancreas/state.yaml" \
  --rename_keys 'input_dataset:output_dataset,input_solution:output_solution' \
  --settings '{"output_scores": "scores.yaml", "output_dataset_info": "dataset_info.yaml", "output_method_configs": "method_configs.yaml", "output_metric_configs": "metric_configs.yaml", "output_task_info": "task_info.yaml", "method_ids": ["bbknn", "mnnpy", "mnnr"]}' \
  --publish_dir "$OUTPUT_DIR" \
  --output_state "state.yaml"

cp trace.txt "$OUTPUT_DIR/trace.txt"


echo ">> Extracting method info"
viash run src/common/process_task_results/get_method_info/config.vsh.yaml -- --input "$OUTPUT_DIR/method_configs.yaml" --output "$OUTPUT_DIR/method_info.json"


echo ">> Uploading results to S3"
aws s3 sync --profile op \
  "resources_test/common/task_metadata/" \
  "s3://openproblems-data/resources_test/common/task_metadata/" \
  --delete --dryrun
