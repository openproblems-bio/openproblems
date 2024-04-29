#!/bin/bash

# Simulating spot-resolution spatial data with alpha = 1

cat > /tmp/params.yaml << 'HERE'
id: spatial_decomposition_process_datasets
input_states: s3://openproblems-data/resources/datasets/**/state.yaml
settings: '{"output_spatial_masked": "$id/spatial_masked.h5ad", "output_single_cell": "$id/single_cell_ref.h5ad", "output_solution": "$id/solution.h5ad", "alpha": 1.0, "simulated_data": "$id/dataset_simulated.h5ad"}'
rename_keys: 'input:output_dataset'
output_state: "$id/state.yaml"
publish_dir: s3://openproblems-data/resources/spatial_decomposition/datasets
HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
  withName:'.*publishStatesProc' {
    memory = '16GB'
    disk = '100GB'
  }
  withLabel:highmem {
    memory = '350GB'
  }
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/spatial_decomposition/workflows/process_datasets/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config /tmp/nextflow.config \
  # --labels spatial_decomposition,process_datasets
