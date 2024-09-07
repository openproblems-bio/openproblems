#!/bin/bash

# only process the 'log_cp10k' datasets
cat > /tmp/params.yaml << 'HERE'
id: predict_modality_process_datasets
input_states: s3://openproblems-data/resources/datasets/**/log_cp10k/state.yaml
settings: '{"output_train_mod1": "$id/train_mod1.h5ad", "output_train_mod2": "$id/train_mod2.h5ad", "output_test_mod1": "$id/test_mod1.h5ad", "output_test_mod2": "$id/test_mod2.h5ad"}'
rename_keys: 'input_mod1:output_mod1,input_mod2:output_mod2'
output_state: "$id/state.yaml"
publish_dir: s3://openproblems-data/resources/predict_modality/datasets
HERE

tw launch https://github.com/openproblems-bio/openproblems.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/predict_modality/workflows/process_datasets/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config src/wf_utils/labels_tw.config \
  --labels predict_modality,process_datasets