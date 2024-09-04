#!/bin/bash

RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"
publish_dir="s3://openproblems-data/resources/denoising/results/${RUN_ID}"

# make sure only log_cp10k is used
cat > /tmp/params.yaml << HERE
input_states: s3://openproblems-data/resources/denoising/datasets/**/log_cp10k/state.yaml
rename_keys: 'input_train:output_train,input_test:output_test'
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/denoising/workflows/run_benchmark/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config src/wf_utils/labels_tw.config \
  --labels denoising,full