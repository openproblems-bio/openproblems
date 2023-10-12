#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

COMMON_DATASETS="resources/datasets/openproblems_v1_multimodal"
OUTPUT_DIR="resources/match_modalities/datasets/openproblems_v1_multimodal"

export NXF_VER=22.04.5

nextflow run . \
  -main-script target/nextflow/match_modalities/workflows/process_datasets/main.nf \
  -profile docker \
  -entry auto \
  -resume \
  --input_states "$COMMON_DATASETS/**/state.yaml" \
  --rename_keys 'input_mod1:output_dataset_mod1,input_mod2:output_dataset_mod2' \
  --settings '{"output_mod1": "$id/output_mod1.h5ad", "output_mod2": "$id/output_mod2.h5ad", "output_solution_mod1": "$id/output_solution_mod1.h5ad", "output_solution_mod2": "$id/output_solution_mod2.h5ad"}' \
  --publish_dir "$OUTPUT_DIR" \
  --output_state '$id/state.yaml'
# output_state should be moved to settings once workaround is solved