#!/bin/bash

cat > /tmp/params.yaml << 'HERE'
input_states: s3://openproblems-data/resources_test/denoising/**/state.yaml
rename_keys: 'input_train:output_train,input_test:output_test'
output_state: "state.yaml"
publish_dir: s3://openproblems-nextflow/temp/denoising/
HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/denoising/workflows/run_benchmark/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config /tmp/nextflow.config \
  --labels denoising,test