#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

OUT_DIR="resources"

echo ">>> Fetching raw results..."
aws s3 sync --profile op \
  s3://openproblems-data/resources/ \
  "$OUT_DIR/" \
  --exclude "*" \
  --include "**/results/run_*/*" \
  --delete

echo ">>> Patch state.yaml files..."
# fix state.yaml id and output_trace
python <<HERE
import os
import re
import glob

def update_state_file(file_path, new_id):
    with open(file_path, 'r') as file:
        content = file.read()
    
    # if output_trace is missing, add it
    if 'output_trace:' not in content:
        content += "\noutput_trace: !file trace.txt\n"
      
    # replace the id with the value of the glob ** pattern
    content = re.sub(r'id: .+', f'id: {new_id}/processed', content)

    with open(file_path, 'w') as file:
        file.write(content)

# find all state.yaml files
state_files = glob.glob('resources/**/state.yaml', recursive=True)
for state_file in state_files:
    # extract the id from the path
    match = re.search(r'resources/(.+?)/state\.yaml', state_file)
    if match:
        new_id = match.group(1)
        update_state_file(state_file, new_id)
        print(f"Updated {state_file} with id: {new_id}")
    else:
        print(f"Could not extract id from {state_file}, skipping.")
HERE

echo ">>> Creating params.yaml..."
cat > /tmp/params.yaml << HERE
input_states: resources/**/state.yaml
rename_keys: 'input_task_info:output_task_info;input_dataset_info:output_dataset_info;input_method_configs:output_method_configs;input_metric_configs:output_metric_configs;input_scores:output_scores;input_trace:output_trace'
output_state: '\$id/state.yaml'
settings: '{"output_combined": "\$id/output_combined.json", "output_report": "\$id/output_report.html", "output_task_info": "\$id/output_task_info.json", "output_dataset_info": "\$id/output_dataset_info.json", "output_method_info": "\$id/output_method_info.json", "output_metric_info": "\$id/output_metric_info.json", "output_results": "\$id/output_results.json", "output_scores": "\$id/output_quality_control.json"}'
publish_dir: "$OUT_DIR"
HERE

echo ">>> Processing results..."
nextflow run target/nextflow/reporting/process_task_results/main.nf \
  -profile docker \
  -params-file /tmp/params.yaml \
  -c common/nextflow_helpers/labels_ci.config \
  -entry auto \
  -resume

# find all files in $OUT with the pattern output_report.html
echo ">>> List reports..."
find "$OUT_DIR" -name "output_report.html"

# echo ">>> Uploading processed results to S3..."
# aws s3 sync --profile op \
#   "resources_test/openproblems/task_results_v4/" \
#   "s3://openproblems-data/resources_test/openproblems/task_results_v4/" \
#   --delete --dryrun

# echo
# echo ">>> Done!"
