#!/bin/bash

DATASET_DIR=resources_test/batch_integration/pancreas

# try running on nf tower
cat > /tmp/params.yaml << 'HERE'
input_states: s3://openproblems-data/resources_test/batch_integration/**/state.yaml
rename_keys: 'input_dataset:output_dataset,input_solution:output_solution'
settings: '{"output": "scores.tsv"}'
output_state: "state.yaml"
publish_dir: s3://openproblems-nextflow/temp/batch_integration/
HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision main_build \
  --pull-latest \
  --main-script src/tasks/batch_integration/workflows/run/main.nf \
  --workspace 53907369739130 \
  --compute-env 1pK56PjjzeraOOC2LDZvN2 \
  --params-file /tmp/params.yaml \
  --config /tmp/nextflow.config