#!/bin/bash
#
#make sure the following command has been executed
#bin/viash_build -q 'label_projection|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

DATASET_DIR=resources_test/common/pancreas

mkdir -p $DATASET_DIR

bin/viash run src/common/dataset_loader/download/config.vsh.yaml -- \
    --url "https://ndownloader.figshare.com/files/24539828" \
    --obs_celltype "celltype" \
    --obs_batch "tech" \
    --name "pancreas" \
    --layer_counts "counts" \
    --output $DATASET_DIR/dataset.h5ad