#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

DATASETS_DIR="resources_test/common"
OUTPUT_DIR="resources_test/predict_modality"

export NXF_VER=22.04.5

echo "Preprocess datasets"
nextflow run . \
  -main-script target/nextflow/predict_modality/workflows/process_datasets/main.nf \
  -profile docker \
  -entry auto \
  -c src/wf_utils/labels_ci.config \
  --input_states "$DATASETS_DIR/**/state.yaml" \
  --rename_keys 'input_rna:output_dataset_rna,input_other_mod:output_dataset_other_mod' \
  --settings '{"output_train_mod1": "$id/train_mod1.h5ad", "output_train_mod2": "$id/train_mod2.h5ad", "output_test_mod1": "$id/test_mod1.h5ad", "output_test_mod2": "$id/test_mod2.h5ad"}' \
  --publish_dir "$OUTPUT_DIR" \
  --output_state '$id/state.yaml'

  echo "Run one method"
  viash run src/tasks/predict_modality/methods/knnr_py/config.vsh.yaml -- \
    --input_train_mod1 $OUTPUT_DIR/bmmc_cite_starter/train_mod1.h5ad \
    --input_train_mod2 $OUTPUT_DIR/bmmc_cite_starter/train_mod2.h5ad \
    --input_test_mod1 $OUTPUT_DIR/bmmc_cite_starter/test_mod1.h5ad \
    --output $OUTPUT_DIR/bmmc_cite_starter/prediction.h5ad
