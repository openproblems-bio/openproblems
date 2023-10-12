#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

COMMON_DATASETS="resources/datasets/openproblems_v1_multimodal"
OUTPUT_DIR="resources/predict_modality/datasets/openproblems_v1_multimodal"

export NXF_VER=22.04.5

nextflow run . \
  -main-script target/nextflow/predict_modality/workflows/process_datasets/main.nf \
  -profile docker \
  -entry auto \
  -resume \
  --input_states "$COMMON_DATASETS/**/state.yaml" \
  --rename_keys 'input_rna:output_dataset_rna,input_other_mod:output_dataset_other_mod' \
  --settings '{"output_train_mod1": "$id/train_mod1.h5ad", "output_train_mod2": "$id/train_mod2.h5ad", "output_test_mod1": "$id/test_mod1.h5ad", "output_test_mod2": "$id/test_mod2.h5ad"}' \
  --publish_dir "$OUTPUT_DIR" \
  --output_state '$id/state.yaml'
# output_state should be moved to settings once workaround is solved