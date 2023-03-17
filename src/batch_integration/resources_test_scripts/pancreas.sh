#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

RAW_DATA=resources_test/common/pancreas/dataset.h5ad
DATASET_DIR=resources_test/batch_integration

if [ ! -f $RAW_DATA ]; then
  echo "Error! Could not find raw data"
  exit 1
fi

mkdir -p $DATASET_DIR

# process dataset
echo process data...
viash run src/batch_integration/split_dataset/config.vsh.yaml -- \
  --input $RAW_DATA \
  --output $DATASET_DIR/unintegrated.h5ad \
  --hvgs 100

echo Running BBKNN
viash run src/batch_integration/methods/bbknn/config.vsh.yaml -- \
  --input $DATASET_DIR/unintegrated.h5ad \
  --output $DATASET_DIR/bbknn.h5ad

echo Running SCVI
viash run src/batch_integration/methods/scvi/config.vsh.yaml -- \
  --input $DATASET_DIR/unintegrated.h5ad \
  --output $DATASET_DIR/scvi.h5ad

echo Running combat
viash run src/batch_integration/methods/combat/config.vsh.yaml -- \
  --input $DATASET_DIR/unintegrated.h5ad \
  --output $DATASET_DIR/combat.h5ad

# run one metric
echo run metrics...
