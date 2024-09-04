#!/bin/bash

# viash run src/common/get_git_sha/config.vsh.yaml -p native -- --input /home/kai/Documents/openroblems/openproblems --output output/op_git_sha.json

TASK_IDS=`ls src/tasks`

for task_id in $TASK_IDS; do
  echo ">> Processing $task_id"
  viash run src/common/get_method_info/config.vsh.yaml -- --input . --task_id $task_id --output output/${task_id}_method.json
  viash run src/migration/check_migration_status/config.vsh.yaml -p native -- --git_sha resources_test/input_git_sha.json --comp_info output/${task_id}_method.json --output output/${task_id}_method_status.json
  viash run src/common/get_metric_info/config.vsh.yaml -- --input . --task_id $task_id --output output/${task_id}_metric.json
  viash run src/migration/check_migration_status/config.vsh.yaml -p native -- --git_sha resources_test/input_git_sha.json --comp_info output/${task_id}_metric.json --output output/${task_id}_metric_status.json

done 