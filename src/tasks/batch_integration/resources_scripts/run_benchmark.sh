#!/bin/bash

RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"
publish_dir="s3://openproblems-data/resources/batch_integration/results/${RUN_ID}"

cat > /tmp/params.yaml << HERE
input_states: s3://openproblems-data/resources/batch_integration/datasets/**/state.yaml
rename_keys: 'input_dataset:output_dataset,input_solution:output_solution'
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/batch_integration/workflows/run_benchmark/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config src/wf_utils/labels_tw.config \
  --labels batch_integration,full