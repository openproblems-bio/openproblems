#!/bin/bash

cat > /tmp/params.yaml << HERE
id: predict_modality_process_datasets
input_states: s3://openproblems-nextflow/resources/datasets/openproblems_v1_multimodal/**/state.yaml
rename_keys: 'input_rna:output_dataset_rna,input_other_mod:output_dataset_other_mod'
settings: '{"output_train_mod1": "train_mod1.h5ad", "output_train_mod2": "train_mod2.h5ad", "output_test_mod1": "test_mod1.h5ad", "output_test_mod2": "test_mod2.h5ad"}'
publish_dir: s3://openproblems-nextflow/resources/predict_modality/datasets/openproblems_v1
HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/predict_modality/workflows/process_datasets/main.nf \
  --workspace 53907369739130 \
  --compute-env 7IkB9ckC81O0dgNemcPJTD \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config /tmp/nextflow.config