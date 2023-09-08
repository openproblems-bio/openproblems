#!/bin/bash

DATASET_DIR=resources_test/common/multimodal

# try running on nf tower
cat > /tmp/params.yaml << HERE
id: scicar
input_mod1: s3://openproblems-data/$DATASET_DIR/dataset_mod1.h5ad
input_mod2: s3://openproblems-data/$DATASET_DIR/dataset_mod2.h5ad
output: scores.tsv
publish_dir: s3://openproblems-nextflow/output_test/v2/match_modalities
HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision integration_build \
  --pull-latest \
  --main-script src/tasks/match_modalities/workflows/run/main.nf \
  --workspace 53907369739130 \
  --compute-env 7IkB9ckC81O0dgNemcPJTD \
  --params-file /tmp/params.yaml \
  --config /tmp/nextflow.config