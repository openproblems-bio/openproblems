#!/bin/bash
#make sure the following command has been executed
#viash ns build -q 'dimensionality_reduction|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

RAW_DATA=resources_test/common/pancreas/dataset.h5ad
DATASET_DIR=resources_test/dimensionality_reduction/pancreas

if [ ! -f $RAW_DATA ]; then
    echo "Error! Could not find raw data"
    exit 1
fi

mkdir -p $DATASET_DIR

# split dataset
viash run src/dimensionality_reduction/split_dataset/config.vsh.yaml -- \
    --input $RAW_DATA \
    --output_train $DATASET_DIR/train.h5ad \
    --output_test $DATASET_DIR/test.h5ad


# run one method
viash run src/dimensionality_reduction/methods/densmap/config.vsh.yaml -- \
    --input $DATASET_DIR/train.h5ad \
    --output $DATASET_DIR/reduced.h5ad

# run one metric
viash run src/dimensionality_reduction/metrics/rmse/config.vsh.yaml -- \
    --input_reduced $DATASET_DIR/reduced.h5ad \
    --input_test $DATASET_DIR/test.h5ad \
    --output $DATASET_DIR/score.h5ad

# run benchmark
export NXF_VER=22.04.5

# after having added a split dataset component
bin/nextflow \
  run . \
  -main-script src/dimensionality_reduction/workflows/run/main.nf \
  -profile docker \
  --id pancreas \
  --dataset_id pancreas \
  --normalization_id log_cpm \
  --input $DATASET_DIR/train.h5ad \
  --input_test $DATASET_DIR/test.h5ad \
  --output scores.tsv \
  --publish_dir $DATASET_DIR/