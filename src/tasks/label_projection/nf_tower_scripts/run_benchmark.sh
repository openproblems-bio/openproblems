#!/bin/bash

cat > /tmp/params.yaml << HERE
id: label_projection
input_states: s3://openproblems-nextflow/resources/label_projection/datasets/**/state.yaml
rename_keys: 'input_train:output_train,input_test:output_test,input_solution:output_solution'
settings: '{"output": "scores.tsv"}'
publish_dir: s3://openproblems-nextflow/output/v2/label_projection
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
  --compute-env 7IkB9ckC81O0dgNemcPJTD \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config /tmp/nextflow.config