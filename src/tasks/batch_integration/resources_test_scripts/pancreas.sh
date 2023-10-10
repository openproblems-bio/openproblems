#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

RAW_DATA=resources_test/common
DATASET_DIR=resources_test/batch_integration

mkdir -p $DATASET_DIR

# process dataset
echo Running process_dataset
nextflow run . \
  -main-script target/nextflow/batch_integration/workflows/process_datasets/main.nf \
  -profile docker \
  -entry auto \
  --input_states "$RAW_DATA/**/state.yaml" \
  --rename_keys 'input:output_dataset' \
  --settings '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad"}' \
  --publish_dir "$DATASET_DIR" \
  --output_state '$id/state.yaml'
# output_state should be moved to settings once workaround is solved

echo Running BBKNN
viash run src/tasks/batch_integration/methods/bbknn/config.vsh.yaml -- \
  --input $DATASET_DIR/pancreas/dataset.h5ad \
  --output $DATASET_DIR/pancreas/integrated_graph.h5ad

echo Running SCVI
viash run src/tasks/batch_integration/methods/scvi/config.vsh.yaml -- \
  --input $DATASET_DIR/pancreas/dataset.h5ad \
  --output $DATASET_DIR/pancreas/integrated_embedding.h5ad

echo Running combat
viash run src/tasks/batch_integration/methods/combat/config.vsh.yaml -- \
  --input $DATASET_DIR/pancreas/dataset.h5ad \
  --output $DATASET_DIR/pancreas/integrated_feature.h5ad