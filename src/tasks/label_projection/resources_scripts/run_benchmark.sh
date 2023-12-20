#!/bin/bash

cat > /tmp/params.yaml << 'HERE'
input_states: s3://openproblems-data/resources/label_projection/datasets/**/state.yaml
rename_keys: 'input_train:output_train,input_test:output_test,input_solution:output_solution'
settings: '{"output": "scores.tsv"}'
output_state: "state.yaml"
publish_dir: s3://openproblems-data/resources/label_projection/results
HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/label_projection/workflows/run_benchmark/main.nf \
  --workspace 53907369739130 \
  --compute-env 1pK56PjjzeraOOC2LDZvN2 \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config /tmp/nextflow.config \
  --labels label_projection,full