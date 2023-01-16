#!/bin/bash



# run a couple of components to generate experimental website view

viash run src/common/list_git_shas/config.vsh.yaml -- \
  --input "../openproblems/" \
  --output "openproblems_git.json"
  
# task_id="label_projection"
# task_id="dimensionality_reduction"
# task_id="denoising"
for task_id in label_projection dimensionality_reduction denoising; do
  out_dir="../website-experimental/results_v2/$task_id/data"

  mkdir -p $out_dir

  # generate task info
  viash run src/common/get_task_info/config.vsh.yaml -- \
    --input "../openproblems-v2" \
    --task_id $task_id \
    --output "$out_dir/task_info.json"

  # generate method info
  viash run src/common/get_method_info/config.vsh.yaml -- \
    --input "../openproblems-v2" \
    --task_id $task_id \
    --output "$out_dir/temp_method_info.json"
  viash run src/common/check_migration_status/config.vsh.yaml -- \
    --git_sha "openproblems_git.json" \
    --comp_info "$out_dir/temp_method_info.json" \
    --output "$out_dir/method_info.json"
  rm "$out_dir/temp_method_info.json"

  # generate metric info
  viash run src/common/get_metric_info/config.vsh.yaml -- \
    --input "../openproblems-v2" \
    --task_id $task_id \
    --output "$out_dir/temp_metric_info.json"
  viash run src/common/check_migration_status/config.vsh.yaml -- \
    --git_sha "openproblems_git.json" \
    --comp_info "$out_dir/temp_metric_info.json" \
    --output "$out_dir/metric_info.json"
  rm "$out_dir/temp_metric_info.json"
  
  # generate results
  viash run src/common/get_results/config.vsh.yaml -- \
    --input_scores "resources/$task_id/benchmarks/openproblems_v1/combined.extract_scores.output.tsv" \
    --input_execution "resources/$task_id/benchmarks/openproblems_v1/nextflow_log.tsv" \
    --output "$out_dir/results.json"
  
  # generate api info
  viash run src/common/get_api_info/config.vsh.yaml -- \
    --input "../openproblems-v2" \
    --task_id $task_id \
    --output "$out_dir/api.json"
done