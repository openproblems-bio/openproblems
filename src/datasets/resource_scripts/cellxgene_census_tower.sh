#!/bin/bash

cat > "/tmp/params.yaml" << 'HERE'
param_list:
  - id: cxg_mouse_pancreas_atlas
    species: mus_musculus
    census_version: "2023-07-25"
    obs_value_filter: "dataset_id == '49e4ffcc-5444-406d-bdee-577127404ba8'"
    obs_batch: donor_id
    dataset_name: Mouse pancreatic islet 
    dataset_summary: Mouse pancreatic islet scRNA-seq atlas across sexes, ages, and stress conditions including diabetes
    dataset_description: To better understand pancreatic β-cell heterogeneity we generated a mouse pancreatic islet atlas capturing a wide range of biological conditions. The atlas contains scRNA-seq datasets of over 300,000 mouse pancreatic islet cells, of which more than 100,000 are β-cells, from nine datasets with 56 samples, including two previously unpublished datasets. The samples vary in sex, age (ranging from embryonic to aged), chemical stress, and disease status (including T1D NOD model development and two T2D models, mSTZ and db/db) together with different diabetes treatments. Additional information about data fields is available in anndata uns field 'field_descriptions' and on https://github.com/theislab/mm_pancreas_atlas_rep/blob/main/resources/cellxgene.md.
    dataset_url: https://cellxgene.cziscience.com/collections/296237e2-393d-4e31-b590-b03f74ac5070
    dataset_reference: hrovatin2023delineating
    dataset_organism: mus_musculus
  - id: cxg_hcla
    species: homo_sapiens
    census_version: "2023-07-25"
    obs_value_filter: "dataset_id == '066943a2-fdac-4b29-b348-40cede398e4e'"
    obs_batch: donor_id
    dataset_name: Human Lung Cell Atlas
    dataset_summary: An integrated cell atlas of the human lung in health and disease (core)
    dataset_description: The integrated Human Lung Cell Atlas (HLCA) represents the first large-scale, integrated single-cell reference atlas of the human lung. It consists of over 2 million cells from the respiratory tract of 486 individuals, and includes 49 different datasets. It is split into the HLCA core, and the extended or full HLCA. The HLCA core includes data of healthy lung tissue from 107 individuals, and includes manual cell type annotations based on consensus across 6 independent experts, as well as demographic, biological and technical metadata.
    dataset_url: https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293
    dataset_reference: sikkema2023integrated
    dataset_organism: homo_sapiens
  - id: cxg_tabula_sapiens
    species: homo_sapiens
    census_version: "2023-07-25"
    obs_value_filter: "dataset_id == '53d208b0-2cfd-4366-9866-c3c6114081bc'"
    obs_batch: [donor_id, assay]
    dataset_name: Tabula Sapiens
    dataset_summary: A multiple-organ, single-cell transcriptomic atlas of humans
    dataset_description: Tabula Sapiens is a benchmark, first-draft human cell atlas of nearly 500,000 cells from 24 organs of 15 normal human subjects. This work is the product of the Tabula Sapiens Consortium. Taking the organs from the same individual controls for genetic background, age, environment, and epigenetic effects and allows detailed analysis and comparison of cell types that are shared between tissues. Our work creates a detailed portrait of cell types as well as their distribution and variation in gene expression across tissues and within the endothelial, epithelial, stromal and immune compartments.
    dataset_url: https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5
    dataset_reference: consortium2022tabula
    dataset_organism: homo_sapiens
  - id: cxg_immune_cell_atlas
    species: homo_sapiens
    census_version: "2023-07-25"
    obs_value_filter: "dataset_id == '1b9d8702-5af8-4142-85ed-020eb06ec4f6'"
    obs_batch: donor_id
    dataset_name: Immune Cell Atlas
    dataset_summary: Cross-tissue immune cell analysis reveals tissue-specific features in humans
    dataset_description: Despite their crucial role in health and disease, our knowledge of immune cells within human tissues remains limited. We surveyed the immune compartment of 16 tissues from 12 adult donors by single-cell RNA sequencing and VDJ sequencing generating a dataset of ~360,000 cells. To systematically resolve immune cell heterogeneity across tissues, we developed CellTypist, a machine learning tool for rapid and precise cell type annotation. Using this approach, combined with detailed curation, we determined the tissue distribution of finely phenotyped immune cell types, revealing hitherto unappreciated tissue-specific features and clonal architecture of T and B cells. Our multitissue approach lays the foundation for identifying highly resolved immune cell types by leveraging a common reference dataset, tissue-integrated expression analysis, and antigen receptor sequencing.
    dataset_url: https://cellxgene.cziscience.com/collections/62ef75e4-cbea-454e-a0ce-998ec40223d3
    dataset_reference: dominguez2022crosstissue
    dataset_organism: homo_sapiens
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
  --revision 62e5f2edeb833e3c932d8ceb89842af79ea3dc1a \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_cellxgene_census/main.nf \
  --workspace 53907369739130 \
  --compute-env 7IkB9ckC81O0dgNemcPJTD \
  --params-file "/tmp/params.yaml" \
  --config /tmp/nextflow.config
