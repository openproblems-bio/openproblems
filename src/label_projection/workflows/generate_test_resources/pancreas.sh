#!/bin/bash
#
#make sure the following command has been executed
#bin/viash_build -q 'label_projection|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

DATASET_DIR=src/label_projection/resources/pancreas

mkdir -p $DATASET_DIR

target/docker/common/data_loader/data_loader\
    --url "https://ndownloader.figshare.com/files/24539828"\
    --name "pancreas"\
    --output $DATASET_DIR/raw_data.h5ad

target/docker/label_projection/data_processing/subsample/subsample\
    --input $DATASET_DIR/raw_data.h5ad\
    --celltype_categories "0:3"\
    --tech_categories "0:-3:-2"\
    --output $DATASET_DIR/toy_data.h5ad

target/docker/label_projection/data_processing/randomize/randomize\
    --input $DATASET_DIR/toy_data.h5ad\
    --output $DATASET_DIR/toy_preprocessed_data.h5ad
