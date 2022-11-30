#!/bin/bash
#
#make sure the following command has been executed
#bin/viash_build -q 'dimensionality_reduction|common'

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
# TODO: implement
# bin/viash run src/dimensionality_reduction/split_dataset/config.vsh.yaml -- \
#     --input $RAW_DATA \
#     --output_dataset $DATASET_DIR/dataset.h5ad \
#     --output_solution $DATASET_DIR/solution.h5ad \
#     --seed 123
cp $RAW_DATA $DATASET_DIR/dataset.h5ad
cp $RAW_DATA $DATASET_DIR/solution.h5ad

# run one method
bin/viash run src/dimensionality_reduction/methods/densmap/config.vsh.yaml -- \
    --input $DATASET_DIR/dataset.h5ad \
    --output $DATASET_DIR/densmap.h5ad

# run one metric
bin/viash run src/dimensionality_reduction/metrics/rmse/config.vsh.yaml -- \
    --input_prediction $DATASET_DIR/densmap.h5ad \
    --input_solution $DATASET_DIR/solution.h5ad \
    --output $DATASET_DIR/densmap_rmse.h5ad

# run benchmark
export NXF_VER=22.04.5

# after having added a split dataset component
# bin/nextflow \
#   run . \
#   -main-script src/dimensionality_reduction/workflows/run/main.nf \
#   -profile docker \
#   -resume \
#   --id pancreas \
#   --dataset_id pancreas \
#   --normalization_id log_cpm \
#   --input_dataset $DATASET_DIR/dataset.h5ad \
#   --input_solution $DATASET_DIR/solution.h5ad \
#   --output scores.tsv \
#   --publish_dir $DATASET_DIR/

bin/nextflow \
  run . \
  -main-script src/dimensionality_reduction/workflows/run/main.nf \
  -profile docker \
  -resume \
  --id pancreas \
  --dataset_id pancreas \
  --normalization_id log_cpm \
  --input $DATASET_DIR/dataset.h5ad \
  --output scores.tsv \
  --publish_dir $DATASET_DIR/