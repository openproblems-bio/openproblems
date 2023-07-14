#!/bin/bash

DATASET_DIR=resources_test/batch_integration/pancreas

# try running on nf tower
cat > /tmp/params.yaml << HERE
id: pancreas_subsample
input: s3://openproblems-data/$DATASET_DIR/unintegrated.h5ad
output: scores.tsv
publish_dir: s3://openproblems-nextflow/output_test/v2/batch_integration
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
  --compute-env 7IkB9ckC81O0dgNemcPJTD \
  --params-file /tmp/params.yaml \
  --config /tmp/nextflow.config