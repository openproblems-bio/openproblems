#!/bin/bash

# Run this prior to executing this script:
# bin/viash_build -q 'spatial_decomposition'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# nextflow run . \
#   -main-script target/nextflow/spatial_decomposition/workflows/process_datasets/main.nf \
#   -profile docker \
#   -entry auto \
#   -c src/wf_utils/labels_ci.config \
#   --id run_test \
#   --input_states "resources_test/common/**/state.yaml" \
#   --rename_keys 'input:output_dataset' \
#   --settings '{"output_spatial_masked": "$id/spatial_masked.h5ad", "output_single_cell": "$id/single_cell_ref.h5ad", "output_solution": "$id/solution.h5ad"}' \
#   --publish_dir "resources_test/spatial_decomposition"

# generate spatial dataset
nextflow run . \
  -main-script target/nextflow/spatial_decomposition/workflows/process_datasets/main.nf \
  -profile docker \
  -entry auto \
  -c src/wf_utils/labels_ci.config \
  --input_states "resources_test/common/**/state.yaml" \
  --rename_keys 'input:output_dataset' \
  --settings '{"output_spatial_masked": "$id/spatial_masked.h5ad", "output_single_cell": "$id/single_cell_ref.h5ad", "output_solution": "$id/solution.h5ad", "generate_dataset": true, "alpha": 1.0, "simulated_data": "$id/dataset_simulated.h5ad"}' \
  --publish_dir "resources_test/spatial_decomposition"