#!/bin/bash

run_date=$(date +%Y%m%d)
publish_dir="s3://openproblems-data/resources/match_modalities/results/${run_date}"

cat > /tmp/params.yaml << HERE
id: match_modalities
input_states: s3://openproblems-data/resources/match_modalities/datasets/**/state.yaml
rename_keys: 'input_mod1:output_mod1,input_mod2:output_mod2,input_solution_mod1:output_solution_mod1,input_solution_mod2:output_solution_mod2'
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
  --main-script target/nextflow/match_modalities/workflows/run_benchmark/main.nf \
  --workspace 53907369739130 \
  --compute-env 1pK56PjjzeraOOC2LDZvN2 \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config /tmp/nextflow.config