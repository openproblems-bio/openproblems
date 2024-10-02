#!/bin/bash

DATASET_DIR=resources_test/common


mkdir -p $DATASET_DIR

wget https://raw.githubusercontent.com/theislab/scib/c993ffd9ccc84ae0b1681928722ed21985fb91d1/scib/resources/g2m_genes_tirosh.txt -O $DATASET_DIR/temp_g2m_genes_tirosh_mm.txt
wget https://raw.githubusercontent.com/theislab/scib/c993ffd9ccc84ae0b1681928722ed21985fb91d1/scib/resources/s_genes_tirosh.txt -O $DATASET_DIR/temp_s_genes_tirosh_mm.txt
KEEP_FEATURES=`cat $DATASET_DIR/temp_g2m_genes_tirosh_mm.txt $DATASET_DIR/temp_s_genes_tirosh_mm.txt | paste -sd ":" -`

cat > "/tmp/params.yaml" << HERE
param_list:
  - id: cxg_immune_cell_atlas
    species: homo_sapiens
    census_version: "2023-07-25"
    obs_value_filter: "dataset_id == '1b9d8702-5af8-4142-85ed-020eb06ec4f6' and donor_id in ['A29', 'A31', 'A35']"
    obs_batch: donor_id
    dataset_name: Immune Cell Atlas
    dataset_summary: Cross-tissue immune cell analysis reveals tissue-specific features in humans
    dataset_description: Despite their crucial role in health and disease, our knowledge of immune cells within human tissues remains limited. We surveyed the immune compartment of 16 tissues from 12 adult donors by single-cell RNA sequencing and VDJ sequencing generating a dataset of ~360,000 cells. To systematically resolve immune cell heterogeneity across tissues, we developed CellTypist, a machine learning tool for rapid and precise cell type annotation. Using this approach, combined with detailed curation, we determined the tissue distribution of finely phenotyped immune cell types, revealing hitherto unappreciated tissue-specific features and clonal architecture of T and B cells. Our multitissue approach lays the foundation for identifying highly resolved immune cell types by leveraging a common reference dataset, tissue-integrated expression analysis, and antigen receptor sequencing.
    dataset_url: https://cellxgene.cziscience.com/collections/62ef75e4-cbea-454e-a0ce-998ec40223d3
    dataset_reference: dominguez2022crosstissue
    dataset_organism: homo_sapiens

normalization_methods: [log_cp10k]
n_obs: 600
n_vars: 1500
output_dataset: '\$id/dataset.h5ad'
output_meta: '\$id/dataset_metadata.yaml'
output_state: '\$id/state.yaml'
output_raw: force_null
output_normalized: force_null
output_pca: force_null
output_hvg: force_null
output_knn: force_null
publish_dir: $DATASET_DIR
do_subsample: true
keep_features: '$KEEP_FEATURES'
HERE

nextflow run . \
  -main-script target/nextflow/datasets/workflows/process_cellxgene_census/main.nf \
  -c src/wf_utils/labels_ci.config \
  -profile docker \
  -params-file "/tmp/params.yaml"

rm -r $DATASET_DIR/temp_*
