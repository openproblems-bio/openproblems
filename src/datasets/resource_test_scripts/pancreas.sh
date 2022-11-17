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

# download dataset
bin/viash run src/datasets/loaders/openproblems_v1/config.vsh.yaml -- \
    --id "pancreas" \
    --obs_celltype "celltype" \
    --obs_batch "tech" \
    --layer_counts "counts" \
    --output $DATASET_DIR/temp_dataset_full.h5ad

# subsample
bin/viash run src/datasets/subsample/config.vsh.yaml -- \
    --input $DATASET_DIR/temp_dataset_full.h5ad \
    --keep_celltype_categories "acinar:beta" \
    --keep_batch_categories "celseq:inDrop4:smarter" \
    --output $DATASET_DIR/temp_dataset0.h5ad \
    --seed 123

# run log cpm normalisation
bin/viash run src/datasets/normalization/log_cpm/config.vsh.yaml -- \
    --input $DATASET_DIR/temp_dataset0.h5ad \
    --output $DATASET_DIR/dataset.h5ad

rm -r $DATASET_DIR/temp_*