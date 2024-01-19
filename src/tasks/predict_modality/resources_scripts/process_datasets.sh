#!/bin/bash

cat > /tmp/params.yaml << 'HERE'
id: predict_modality_process_datasets
input_states: s3://openproblems-data/resources/datasets/**/state.yaml
rename_keys: 'input_rna:output_rna,input_other_mod:output_other_mod'
settings: '{"output_train_mod1": "$id/train_mod1.h5ad", "output_train_mod2": "$id/train_mod2.h5ad", "output_test_mod1": "$id/test_mod1.h5ad", "output_test_mod2": "$id/test_mod2.h5ad"}'
output_state: "$id/state.yaml"
publish_dir: s3://openproblems-data/resources/predict_modality/datasets
HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
  withName:'.*publishStatesProc' {
    memory = '16GB'
    disk = '100GB'
  }
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/predict_modality/workflows/process_datasets/main.nf \
  --workspace 53907369739130 \
  --compute-env 1pK56PjjzeraOOC2LDZvN2 \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config /tmp/nextflow.config \
  --labels predict_modality,process_datasets