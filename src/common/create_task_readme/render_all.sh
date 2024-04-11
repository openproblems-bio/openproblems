#!/bin/bash

set -e

TASK_IDS=`ls src/tasks`

for task_id in $TASK_IDS; do
  echo ">> Processing $task_id"
  viash run src/common/create_task_readme/config.vsh.yaml -- --task $task_id
done