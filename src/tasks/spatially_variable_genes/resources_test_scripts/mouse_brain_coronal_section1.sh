#!/bin/bash

# make sure the following command has been executed
# viash ns build -q 'spatially_variable_genes|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

RAW_DATA=resources_test/common
DATASET_DIR=resources_test/spatially_variable_genes

mkdir -p $DATASET_DIR

echo "Running process_dataset"
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

echo "Running control method"
viash run src/tasks/spatially_variable_genes/control_methods/true_ranking/config.vsh.yaml -- \
  --input_data $DATASET_DIR/mouse_brain_coronal_section1/dataset.h5ad \
  --input_solution $DATASET_DIR/mouse_brain_coronal_section1/solution.h5ad \
  --output $DATASET_DIR/mouse_brain_coronal_section1/output.h5ad

echo "Running metric"
viash run src/tasks/spatially_variable_genes/metrics/correlation/config.vsh.yaml -- \
    --input_method $DATASET_DIR/mouse_brain_coronal_section1/output.h5ad \
    --input_solution $DATASET_DIR/mouse_brain_coronal_section1/solution.h5ad \
    --output $DATASET_DIR/mouse_brain_coronal_section1/score.h5ad
