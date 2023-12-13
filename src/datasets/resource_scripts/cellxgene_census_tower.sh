#!/bin/bash

cat > "/tmp/params.yaml" << 'HERE'
param_list:

  # template for adding new datasets
  # - id: cxg_
  #   species:
  #   census_version: "2023-07-25"
  #   obs_value_filter: "dataset_id == ''"
  #   obs_batch:
  #   dataset_name:
  #   dataset_summary:
  #   dataset_description:
  #   dataset_url:
  #   dataset_reference:
  #   dataset_organism:

normalization_methods: [log_cp10k, sqrt_cp10k, l1_sqrt]
output_dataset: '$id/dataset.h5ad'
output_meta: '$id/dataset_metadata.yaml'
output_state: '$id/state.yaml'
output_raw: force_null
output_normalized: force_null
output_pca: force_null
output_hvg: force_null
output_knn: force_null
publish_dir: s3://openproblems-nextflow/resources/datasets/cellxgene_census
HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_cellxgene_census/main.nf \
  --workspace 53907369739130 \
  --compute-env 7IkB9ckC81O0dgNemcPJTD \
  --params-file "/tmp/params.yaml" \
  --config /tmp/nextflow.config
