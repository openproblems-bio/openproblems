#!/bin/bash

RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"
publish_dir="s3://openproblems-data/resources/spatial_decomposition/results/${RUN_ID}"

cat > /tmp/params.yaml << HERE
input_states: s3://openproblems-data/resources/spatial_decomposition/datasets/**/state.yaml
rename_keys: 'input_single_cell:output_single_cell,input_spatial_masked:output_spatial_masked,input_solution:output_solution'
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/spatial_decomposition/workflows/run_benchmark/main.nf \
  --workspace 53907369739130 \
  --compute-env 1pK56PjjzeraOOC2LDZvN2 \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config src/wf_utils/labels_tw.config \
  --labels spatial_decomposition,full