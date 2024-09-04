#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

RAW_DATA=resources_test/common
DATASET_DIR=resources_test/label_projection

mkdir -p $DATASET_DIR

# process dataset
echo Running process_dataset
nextflow run . \
  -main-script target/nextflow/label_projection/workflows/process_datasets/main.nf \
  -profile docker \
  -entry auto \
  --input_states "$RAW_DATA/**/state.yaml" \
  --rename_keys 'input:output_dataset' \
  --settings '{"output_train": "$id/train.h5ad", "output_test": "$id/test.h5ad", "output_solution": "$id/solution.h5ad"}' \
  --publish_dir "$DATASET_DIR" \
  --output_state '$id/state.yaml'
# output_state should be moved to settings once workaround is solved

# run one method
viash run src/tasks/label_projection/methods/knn/config.vsh.yaml -- \
    --input_train $DATASET_DIR/pancreas/train.h5ad \
    --input_test $DATASET_DIR/pancreas/test.h5ad \
    --output $DATASET_DIR/pancreas/prediction.h5ad

# run one metric
viash run src/tasks/label_projection/metrics/accuracy/config.vsh.yaml -- \
    --input_prediction $DATASET_DIR/pancreas/prediction.h5ad \
    --input_solution $DATASET_DIR/pancreas/solution.h5ad \
    --output $DATASET_DIR/pancreas/score.h5ad
