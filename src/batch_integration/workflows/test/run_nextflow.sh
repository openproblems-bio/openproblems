#!/bin/bash

# Run this prior to executing this script:
# bin/viash_build -q 'batch_integration'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -xe

RAW_DATA=resources_test/common/pancreas/dataset.h5ad
DATASET_DIR=resources_test/batch_integration/pancreas

if [ ! -d "$DATASET_DIR" ]; then
  mkdir -p "$DATASET_DIR"
fi

# viash run src/batch_integration/datasets/preprocessing/config.vsh.yaml -- \
#   --input $RAW_DATA \
#   --output $DATASET_DIR/processed.h5ad \
#   --label celltype \
#   --batch tech \
#   --hvgs 100

# choose a particular version of nextflow
# export NXF_VER=21.10.6

# bin/nextflow \
#   run . \
#   -main-script src/batch_integration/workflows/run/main.nf \
#   --publishDir output/batch_integration \
#   -resume \
#   -with-docker


# run benchmark
export NXF_VER=22.04.5

  # -profile docker \
nextflow \
  run . \
  -main-script src/batch_integration/workflows/test/main.nf \
  -with-docker \
  -resume \
  --id pancreas \
  --dataset_id pancreas \
  --input $DATASET_DIR/processed.h5ad \
  --output scores.tsv \
  --publish_dir $DATASET_DIR/
