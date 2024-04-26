#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

RAW_DATA=resources_test/common
DATASET_DIR=resources_test/match_modalities

mkdir -p $DATASET_DIR

# process dataset
echo Running process_dataset
nextflow run . \
  -main-script target/nextflow/match_modalities/workflows/process_datasets/main.nf \
  -profile docker \
  -entry auto \
  --input_states "$RAW_DATA/**/state.yaml" \
  --rename_keys 'input_mod1:output_mod1,input_mod2:output_mod2' \
  --settings '{"output_mod1": "$id/dataset_mod1.h5ad", "output_mod2": "$id/dataset_mod2.h5ad", "output_solution_mod1": "$id/solution_mod1.h5ad", "output_solution_mod2": "$id/solution_mod2.h5ad"}' \
  --publish_dir "$DATASET_DIR" \
  --output_state '$id/state.yaml'
# output_state should be moved to settings once workaround is solved

# run one method
viash run src/tasks/match_modalities/methods/fastmnn/config.vsh.yaml -- \
    --input_mod1 $DATASET_DIR/scicar_cell_lines/dataset_mod1.h5ad \
    --input_mod2 $DATASET_DIR/scicar_cell_lines/dataset_mod2.h5ad \
    --output_mod1 $DATASET_DIR/scicar_cell_lines/integrated_mod1.h5ad \
    --output_mod2 $DATASET_DIR/scicar_cell_lines/integrated_mod2.h5ad
