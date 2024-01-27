#!/bin/bash

set -e

cat > "/tmp/params.yaml" << 'HERE'
param_list:
  - id: cxg_mouse_pancreas_atlas
    species: mus_musculus
    census_version: "2023-07-25"
    obs_value_filter: "dataset_id == '49e4ffcc-5444-406d-bdee-577127404ba8' and donor_id in ['mouse_pancreatic_islet_atlas_Hrovatin__Fltp_2y__MUC13974', 'mouse_pancreatic_islet_atlas_Hrovatin__Fltp_2y__MUC13975', 'mouse_pancreatic_islet_atlas_Hrovatin__Fltp_2y__MUC13976']"
    obs_batch: donor_id
    dataset_name: Mouse Pancreatic Islet Atlas
    dataset_summary: Mouse pancreatic islet scRNA-seq atlas across sexes, ages, and stress conditions including diabetes
    dataset_description: To better understand pancreatic β-cell heterogeneity we generated a mouse pancreatic islet atlas capturing a wide range of biological conditions. The atlas contains scRNA-seq datasets of over 300,000 mouse pancreatic islet cells, of which more than 100,000 are β-cells, from nine datasets with 56 samples, including two previously unpublished datasets. The samples vary in sex, age (ranging from embryonic to aged), chemical stress, and disease status (including T1D NOD model development and two T2D models, mSTZ and db/db) together with different diabetes treatments. Additional information about data fields is available in anndata uns field 'field_descriptions' and on https://github.com/theislab/mm_pancreas_atlas_rep/blob/main/resources/cellxgene.md.
    dataset_url: https://cellxgene.cziscience.com/collections/296237e2-393d-4e31-b590-b03f74ac5070
    dataset_reference: hrovatin2023delineating
    dataset_organism: mus_musculus

normalization_methods: [log_cp10k]
output_dataset: '$id/dataset.h5ad'
output_meta: '$id/dataset_metadata.yaml'
output_state: '$id/state.yaml'
output_raw: force_null
output_normalized: force_null
output_pca: force_null
output_hvg: force_null
output_knn: force_null
publish_dir: resources_test/common
do_subsample: true
HERE

nextflow run . \
  -main-script target/nextflow/datasets/workflows/process_cellxgene_census/main.nf \
  -profile docker \
  -params-file "/tmp/params.yaml"

# src/tasks/batch_integration/resources_test_scripts/process.sh