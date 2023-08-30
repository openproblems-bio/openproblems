#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

RAW_DATA=resources_test/common/pancreas/dataset.h5ad
DATASET_DIR=resources_test/batch_integration/pancreas

if [ ! -f $RAW_DATA ]; then
  echo "Error! Could not find raw data"
  exit 1
fi

mkdir -p $DATASET_DIR

# process dataset
echo process data...
viash run src/tasks/batch_integration/process_dataset/config.vsh.yaml -- \
  --input $RAW_DATA \
  --output_dataset $DATASET_DIR/dataset.h5ad \
  --output_solution $DATASET_DIR/solution.h5ad \
  --hvgs 100

echo Running BBKNN
viash run src/tasks/batch_integration/methods/bbknn/config.vsh.yaml -- \
  --input $DATASET_DIR/dataset.h5ad \
  --output $DATASET_DIR/integrated_graph.h5ad

echo Running SCVI
viash run src/tasks/batch_integration/methods/scvi/config.vsh.yaml -- \
  --input $DATASET_DIR/dataset.h5ad \
  --output $DATASET_DIR/integrated_embedding.h5ad

echo Running combat
viash run src/tasks/batch_integration/methods/combat/config.vsh.yaml -- \
  --input $DATASET_DIR/dataset.h5ad \
  --output $DATASET_DIR/integrated_feature.h5ad

# run one metric
echo run metrics...
