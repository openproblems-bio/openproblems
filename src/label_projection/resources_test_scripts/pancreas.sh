#!/bin/bash
#
#make sure the following command has been executed
#bin/viash_build -q 'label_projection|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

RAW_DATA=resources_test/common/pancreas/dataset_cpm.h5ad
DATASET_DIR=resources_test/label_projection/pancreas

if [ ! -f $RAW_DATA ]; then
    echo "Error! Could not find raw data"
    exit 1
fi

mkdir -p $DATASET_DIR

# censor dataset
bin/viash run src/label_projection/data_processing/censoring/config.vsh.yaml -- \
    --input $RAW_DATA \
    --output_train $DATASET_DIR/dataset_cpm_train.h5ad \
    --output_test $DATASET_DIR/dataset_cpm_test.h5ad \
    --output_solution $DATASET_DIR/dataset_cpm_solution.h5ad \
    --seed 123

# run one method
bin/viash run src/label_projection/methods/knn_classifier/config.vsh.yaml -- \
    --input_train $DATASET_DIR/dataset_cpm_train.h5ad \
    --input_test $DATASET_DIR/dataset_cpm_test.h5ad \
    --output $DATASET_DIR/dataset_cpm_knn.h5ad

# run one metric
bin/viash run src/label_projection/metrics/accuracy/config.vsh.yaml -- \
    --input_prediction $DATASET_DIR/dataset_cpm_knn.h5ad \
    --input_solution $DATASET_DIR/dataset_cpm_solution.h5ad \
    --output $DATASET_DIR/dataset_cpm_knn_accuracy.h5ad
