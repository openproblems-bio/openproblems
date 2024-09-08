#!/bin/bash

# Run this prior to executing this script:
# viash ns build -q 'spatially_variable_genes'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

RAW_DATA=resources_test/common
DATASET_DIR=resources_test/spatially_variable_genes

mkdir -p $DATASET_DIR

nextflow run . \
  -main-script target/nextflow/spatially_variable_genes/workflows/process_datasets/main.nf \
  -profile docker \
  -c src/wf_utils/labels_ci.config \
  --id mouse_brain_coronal_section1 \
  --input $RAW_DATA/mouse_brain_coronal_section1/dataset.h5ad \
  --output_dataset dataset.h5ad \
  --output_solution solution.h5ad \
  --dataset_simulated_normalized simulated_dataset.h5ad \
  --publish_dir $DATASET_DIR/mouse_brain_coronal_section1 \
  --output_state "state.yaml" \
  --gp_k_sim 50 \
  --select_top_variable_genes 50 \
  --num_reference_genes 200