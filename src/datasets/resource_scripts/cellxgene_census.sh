#!/bin/bash

cat > "/tmp/params.yaml" << 'HERE'
param_list:
  - id: cxg_mm_pancreas_atlas
    species: mus_musculus
    census_version: "2023-07-25"
    obs_value_filter: "dataset_id == '49e4ffcc-5444-406d-bdee-577127404ba8'"
    obs_batch: donor_id
    dataset_name: Mouse pancreatic islet atlas
    dataset_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211799
    dataset_reference: hrovatin2023delineating
    dataset_summary: Mouse pancreatic islet scRNA-seq atlas across sexes, ages, and stress conditions including diabetes
    dataset_description: To better understand pancreatic β-cell heterogeneity we generated a mouse pancreatic islet atlas capturing a wide range of biological conditions. The atlas contains scRNA-seq datasets of over 300,000 mouse pancreatic islet cells, of which more than 100,000 are β-cells, from nine datasets with 56 samples, including two previously unpublished datasets. The samples vary in sex, age (ranging from embryonic to aged), chemical stress, and disease status (including T1D NOD model development and two T2D models, mSTZ and db/db) together with different diabetes treatments. Additional information about data fields is available in anndata uns field 'field_descriptions' and on https://github.com/theislab/mm_pancreas_atlas_rep/blob/main/resources/cellxgene.md.
    dataset_organism: mus_musculus

normalization_methods: [log_cp10k, sqrt_cp10k, l1_sqrt]
output_dataset: '$id/dataset.h5ad'
output_meta: '$id/dataset_metadata.yaml'
output_state: '$id/state.yaml'
output_raw: force_null
output_normalized: force_null
output_pca: force_null
output_hvg: force_null
output_knn: force_null
publish_dir: output/temp
HERE

nextflow run . \
  -main-script target/nextflow/datasets/workflows/process_cellxgene_census/main.nf \
  -profile docker \
  -params-file "/tmp/params.yaml"