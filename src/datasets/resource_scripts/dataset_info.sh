#!/bin/bash

DATASETS_DIR="s3://openproblems-data/resources/datasets"

cat > "/tmp/params.yaml" << HERE
param_list:
  - id: openproblems_v1
    input_states: "$DATASETS_DIR/openproblems_v1/**/log_cp10k/state.yaml"
    rename_keys: 'input:output_dataset'
  - id: openproblems_v1_multimodal
    input_states: "$DATASETS_DIR/openproblems_v1_multimodal/**/log_cp10k/state.yaml"
    rename_keys: 'input:output_mod1'
  - id: cellxgene_census
    input_states: "$DATASETS_DIR/cellxgene_census/**/log_cp10k/state.yaml"
    rename_keys: 'input:output_dataset'
settings: '{"output": "dataset_info.yaml"}'
output_state: state.yaml
publish_dir: "$DATASETS_DIR"
HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
  withLabel: highmem {
    memory = '350GB'
  }
  withName: '.*publishStatesProc' {
    memory = '16GB'
    disk = '100GB'
  }
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision main_build \
  --entry-name auto \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/extract_dataset_info/main.nf \
  --workspace 53907369739130 \
  --compute-env 1pK56PjjzeraOOC2LDZvN2 \
  --params-file "/tmp/params.yaml" \
  --config /tmp/nextflow.config


# # run locally after the above has finished
# nextflow run . \
#   -main-script target/nextflow/common/process_task_results/get_dataset_info/main.nf \
#   -profile docker \
#   -resume \
#   --input "$DATASETS_DIR/dataset_info.yaml" \
#   --task_id "common" \
#   --output "dataset_info.json" \
#   --output_state state.yaml \
#   --publish_dir "../website/documentation/reference/datasets/data/"