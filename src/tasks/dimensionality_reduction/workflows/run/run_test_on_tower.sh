#!/bin/bash

DATASET_DIR=resources_test/dimensionality_reduction/pancreas

# try running on nf tower
cat > /tmp/params.yaml << HERE
id: pancreas_subsample
input: s3://openproblems-data/$DATASET_DIR/dataset.h5ad
input_solution: s3://openproblems-data/$DATASET_DIR/solution.h5ad
dataset_id: pancreas
normalization_id: log_cpm
output: scores.tsv
publish_dir: s3://openproblems-nextflow/output_test/v2/dimensionality_reduction
HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision integration_build \
  --pull-latest \
  --main-script src/tasks/dimensionality_reduction/workflows/run/main.nf \
  --workspace 53907369739130 \
  --compute-env 7IkB9ckC81O0dgNemcPJTD \
  --params-file /tmp/params.yaml \
  --config /tmp/nextflow.config