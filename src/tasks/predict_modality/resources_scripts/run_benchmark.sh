#!/bin/bash

run_date=$(date +%Y%m%d)
publish_dir="s3://openproblems-data/resources/predict_modality/results/${run_date}"

cat > /tmp/params.yaml << HERE
id: predict_modality
input_states: s3://openproblems-data/resources/predict_modality/datasets/**/state.yaml
rename_keys: 'input_train_mod1:output_train_mod1,input_train_mod2:output_train_mod2,input_test_mod1:output_test_mod1,input_test_mod2:output_test_mod2'
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
}

trace {
    enabled = true
    overwrite = true
    file    = "$publish_dir/trace.txt"
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/predict_modality/workflows/run_benchmark/main.nf \
  --workspace 53907369739130 \
  --compute-env 1pK56PjjzeraOOC2LDZvN2 \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config /tmp/nextflow.config