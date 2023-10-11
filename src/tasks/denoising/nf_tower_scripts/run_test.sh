#!/bin/bash

DATASET_DIR=resources_test/denoising/pancreas

# try running on nf tower
cat > /tmp/params.yaml << HERE
id: denoising_test
input_states: s3://openproblems-data/resources_test/denoising/pancreas/
rename_keys: 'input_train:output_train,input_test:output_test'
settings: '{"output": "scores.tsv"}'
publish_dir: s3://openproblems-nextflow/output_test/v2/denoising/
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
  --compute-env 7IkB9ckC81O0dgNemcPJTD \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config /tmp/nextflow.config