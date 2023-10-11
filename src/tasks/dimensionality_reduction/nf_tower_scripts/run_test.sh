#!/bin/bash


# try running on nf tower
cat > /tmp/params.yaml << HERE
id: dimensionality_reduction
input_states: s3://openproblems-data/resources_test/dimensionality_reduction/pancreas
rename_keys: 'input_dataset:output_dataset,input_solution:output_solution'
settings: '{"output": "scores.tsv"}'
publish_dir: s3://openproblems-nextflow/output_test/v2/dimensionality_reduction
HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/dimensionality_reduction/workflows/run_benchmark/main.nf \
  --workspace 53907369739130 \
  --compute-env 7IkB9ckC81O0dgNemcPJTD \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config /tmp/nextflow.config