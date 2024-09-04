#!/bin/bash

set -e

cat > /tmp/params.yaml << 'HERE'
param_list:
  - id: mouse_cerebellum
    input_data: "https://zenodo.org/records/12785822/files/Slide-seqV2_stickels2020highly_stickels2021highly_SlideSeqV2_Mouse_Olfactory_bulb_Puck_200127_15_data_whole.h5ad?download=1"
    dataset_name: Mouse cerebellum
    dataset_url: "..."
    dataset_summary: ...
    dataset_description: "..."
    dataset_reference: ref
    dataset_organism: Mus musculus

normalization_methods: [log_cp10k]
n_obs: 600
n_vars: 500
output_dataset: '$id/dataset.h5ad'
output_meta: '$id/dataset_metadata.yaml'
output_state: '$id/state.yaml'
output_raw: force_null
output_normalized: force_null
publish_dir: resources_test/common
do_subsample: true
spot_filter_min_genes: 200
gene_filter_min_spots: 50
remove_mitochondrial: true
HERE

nextflow run . \
  -main-script target/nextflow/datasets/workflows/process_spatial_from_zenodo/main.nf \
  -c src/wf_utils/labels_ci.config \
  -profile docker \
  -params-file "/tmp/params.yaml"

