#!/bin/bash
#make sure the following command has been executed
#viash ns build -q 'dimensionality_reduction|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

RAW_DATA=resources_test/common
DATASET_DIR=resources_test/dimensionality_reduction

mkdir -p $DATASET_DIR

# process dataset
echo Running process_dataset
nextflow run . \
  -main-script target/nextflow/dimensionality_reduction/workflows/process_datasets/main.nf \
  -profile docker \
  -entry auto \
  --input_states "$RAW_DATA/**/state.yaml" \
  --rename_keys 'input:output_dataset' \
  --settings '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad"}' \
  --publish_dir "$DATASET_DIR" \
  --output_state '$id/state.yaml'
# output_state should be moved to settings once workaround is solved


# run one method
viash run src/tasks/dimensionality_reduction/methods/densmap/config.vsh.yaml -- \
    --input $DATASET_DIR/dataset.h5ad \
    --output $DATASET_DIR/embedding.h5ad

# run one metric
viash run src/tasks/dimensionality_reduction/metrics/distance_correlation/config.vsh.yaml -- \
    --input_embedding $DATASET_DIR/embedding.h5ad \
    --input_solution $DATASET_DIR/solution.h5ad \
    --output $DATASET_DIR/score.h5ad

# # run benchmark
# export NXF_VER=22.04.5

# # after having added a split dataset component
# nextflow \
#   run . \
#   -main-script src/tasks/dimensionality_reduction/workflows/run/main.nf \
#   -profile docker \
#   --id pancreas \
#   --input_dataset $DATASET_DIR/dataset.h5ad \
#   --input_solution $DATASET_DIR/solution.h5ad \
#   --output scores.tsv \
#   --publish_dir $DATASET_DIR/