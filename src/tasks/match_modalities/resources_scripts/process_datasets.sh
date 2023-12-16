#!/bin/bash

cat > /tmp/params.yaml << 'HERE'
id: match_modalities_process_datasets
input_states: s3://openproblems-data/resources/datasets/openproblems_v1_multimodal/**/state.yaml
rename_keys: 'input_mod1:output_dataset_mod1,input_mod2:output_dataset_mod2'
settings: '{"output_mod1": "$id/output_mod1.h5ad", "output_mod2": "$id/output_mod2.h5ad", "output_solution_mod1": "$id/output_solution_mod1.h5ad", "output_solution_mod2": "$id/output_solution_mod2.h5ad"}'
output_state: "$id/state.yaml"
publish_dir: s3://openproblems-data/resources/match_modalities/datasets/openproblems_v1_multimodal
HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/match_modalities/workflows/process_datasets/main.nf \
  --workspace 53907369739130 \
  --compute-env 1pK56PjjzeraOOC2LDZvN2 \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config /tmp/nextflow.config
