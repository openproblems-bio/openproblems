#!/bin/bash

# task_id="label_projection"
# task_id="dimensionality_reduction"
# task_id="denoising"

# run a couple of components to generate experimental website view

for task_id in label_projection dimensionality_reduction denoising; do
  out_dir="../website-experimental/results_v2/$task_id/data"

  mkdir -p $out_dir

  viash run src/common/get_method_info/config.vsh.yaml -p native -- \
    --input "src/$task_id" \
    --output "$out_dir/method_info.json"
  viash run src/common/get_metric_info/config.vsh.yaml -p native -- \
    --input "src/$task_id" \
    --output "$out_dir/metric_info.json"
  viash run src/common/get_results/config.vsh.yaml -p native -- \
    --input_scores "resources/$task_id/benchmarks/openproblems_v1/combined.extract_scores.output.tsv" \
    --input_execution "resources/$task_id/benchmarks/openproblems_v1/nextflow_log.tsv" \
    --output "$out_dir/results.json"
  viash run src/common/get_api_info/config.vsh.yaml -p native -- \
    --input "src/$task_id" \
    --output "$out_dir/api.json"
done