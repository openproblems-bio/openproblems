#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

RAW_DATA=resources_test/common
DATASET_DIR=resources_test/batch_integration

mkdir -p $DATASET_DIR

# process dataset
echo Running process_dataset
nextflow run . \
  -main-script src/tasks/batch_integration/workflows/process_datasets/main.nf \
  -profile docker \
  -entry auto \
  --id resources_test \
  --input_states "$RAW_DATA/**/state.yaml" \
  --rename_keys 'input:output_dataset' \
  --settings '{"output_dataset": "dataset.h5ad", "output_solution": "solution.h5ad"}' \
  --publish_dir "$DATASET_DIR"

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