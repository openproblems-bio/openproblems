#!/bin/bash

cat > /tmp/params.yaml << 'HERE'
id: dimensionality_reduction_process_datasets
input_states: s3://openproblems-data/resources/datasets/**/state.yaml
rename_keys: 'input:output_dataset'
settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad"}'
output_state: "$id/state.yaml"
publish_dir: s3://openproblems-data/resources/dimensionality_reduction/datasets
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

tw launch https://github.com/openproblems-bio/openproblems.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/dimensionality_reduction/workflows/process_datasets/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config /tmp/nextflow.config \
  --labels dimensionality_reduction,process_datasets