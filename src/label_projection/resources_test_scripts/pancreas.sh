#!/bin/bash
#
#make sure the following command has been executed
#viash_build -q 'label_projection|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

RAW_DATA=resources_test/common/pancreas/dataset.h5ad
DATASET_DIR=resources_test/label_projection/pancreas

if [ ! -f $RAW_DATA ]; then
    echo "Error! Could not find raw data"
    exit 1
fi

mkdir -p $DATASET_DIR

# split dataset
viash run src/label_projection/split_dataset/config.vsh.yaml -- \
    --input $RAW_DATA \
    --output_train $DATASET_DIR/train.h5ad \
    --output_test $DATASET_DIR/test.h5ad \
    --output_solution $DATASET_DIR/solution.h5ad \
    --seed 123

# run one method
viash run src/label_projection/methods/knn/config.vsh.yaml -- \
    --input_train $DATASET_DIR/train.h5ad \
    --input_test $DATASET_DIR/test.h5ad \
    --output $DATASET_DIR/knn.h5ad

# run one metric
viash run src/label_projection/metrics/accuracy/config.vsh.yaml -- \
    --input_prediction $DATASET_DIR/knn.h5ad \
    --input_solution $DATASET_DIR/solution.h5ad \
    --output $DATASET_DIR/knn_accuracy.h5ad

# run benchmark
export NXF_VER=22.04.5

nextflow \
  run . \
  -main-script src/label_projection/workflows/run/main.nf \
  -profile docker \
  -resume \
  --id pancreas \
  --dataset_id pancreas \
  --normalization_id log_cpm \
  --input_train $DATASET_DIR/train.h5ad \
  --input_test $DATASET_DIR/test.h5ad \
  --input_solution $DATASET_DIR/solution.h5ad \
  --output scores.tsv \
  --publish_dir $DATASET_DIR/