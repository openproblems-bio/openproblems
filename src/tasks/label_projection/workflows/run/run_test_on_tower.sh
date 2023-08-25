#!/bin/bash

DATASET_DIR=resources_test/label_projection/pancreas

# try running on nf tower
cat > /tmp/params.yaml << HERE
id: pancreas_subsample
input_train: s3://openproblems-data/$DATASET_DIR/train.h5ad
input_test: s3://openproblems-data/$DATASET_DIR/test.h5ad
input_solution: s3://openproblems-data/$DATASET_DIR/solution.h5ad
dataset_id: pancreas
normalization_id: log_cp10k
output: scores.tsv
publish_dir: s3://openproblems-nextflow/output_test/v2/label_projection
HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision main_build \
  --pull-latest \
  --main-script src/tasks/label_projection/workflows/run/main.nf \
  --workspace 53907369739130 \
  --compute-env 7IkB9ckC81O0dgNemcPJTD \
  --params-file /tmp/params.yaml \
  --config /tmp/nextflow.config