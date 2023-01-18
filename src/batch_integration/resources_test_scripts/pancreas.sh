#!/bin/bash
set -xe
#
#make sure the following command has been executed
#bin/viash_build -q 'batch_integration|common'

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

# build components
bin/viash_build -q batch

# process dataset
echo process data...
bin/viash run src/batch_integration/datasets/preprocessing/config.vsh.yaml -- \
  --input $RAW_DATA \
  --output $DATASET_DIR/pancreas/processed.h5ad \
  --label celltype \
  --batch tech \
  --hvgs 100

# run methods
echo run methods...

# Graph methods
bin/viash run src/batch_integration/graph/methods/bbknn/config.vsh.yaml -- \
  --input $DATASET_DIR/pancreas/processed.h5ad \
  --output $DATASET_DIR/graph/methods/bbknn.h5ad

bin/viash run src/batch_integration/graph/methods/combat/config.vsh.yaml -- \
  --input $DATASET_DIR/pancreas/processed.h5ad \
  --output $DATASET_DIR/graph/methods/combat.h5ad

# Embedding method
bin/viash run src/batch_integration/embedding/methods/combat/config.vsh.yaml -- \
  --input $DATASET_DIR/pancreas/processed.h5ad \
  --output $DATASET_DIR/embedding/methods/combat.h5ad

# run one metric
echo run metrics...
