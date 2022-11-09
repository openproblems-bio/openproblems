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
bin/viash run src/common/dataset_loader/download/config.vsh.yaml -- \
    --url "https://ndownloader.figshare.com/files/24539828" \
    --obs_celltype "celltype" \
    --obs_batch "tech" \
    --name "pancreas" \
    --layer_counts "counts" \
    --output $DATASET_DIR/temp_full_dataset.h5ad

# subsample
bin/viash run src/common/subsample/config.vsh.yaml -- \
    --input $DATASET_DIR/temp_full_dataset.h5ad \
    --keep_celltype_categories "acinar:beta" \
    --keep_batch_categories "celseq:inDrop4:smarter" \
    --output $DATASET_DIR/dataset.h5ad \
    --seed 123

# run one normalisation
bin/viash run src/common/normalization/log_cpm/config.vsh.yaml -- \
    --input $DATASET_DIR/dataset.h5ad \
    --output $DATASET_DIR/dataset_cpm.h5ad
