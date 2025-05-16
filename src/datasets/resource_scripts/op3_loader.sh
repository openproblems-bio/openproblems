#!/bin/bash

set -e

params_file="/tmp/datasets_openproblems_neurips2022_params.yaml"

## TODO: fill in the parameters below
echo "TODO: fill in the parameters below" && exit 1

cat > "$params_file" << 'HERE'
param_list:
  - id: openproblems_neurips2022/pbmc_cite
    input: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279945/suppl/GSE279945_sc_counts_processed.h5ad
    dataset_name: OpenProblems NeurIPS2022 CITE-Seq
    dataset_organism: homo_sapiens
    dataset_summary: Single-cell CITE-Seq (GEX+ADT) data collected from bone marrow mononuclear cells of 12 healthy human donors.
    dataset_description: "Single-cell CITE-Seq data collected from bone marrow mononuclear cells of 12 healthy human donors using the 10X 3 prime Single-Cell Gene Expression kit with Feature Barcoding in combination with the BioLegend TotalSeq B Universal Human Panel v1.0. The dataset was generated to support Multimodal Single-Cell Data Integration Challenge at NeurIPS 2022. Samples were prepared using a standard protocol at four sites. The resulting data was then annotated to identify cell types and remove doublets. The dataset was designed with a nested batch layout such that some donor samples were measured at multiple sites with some donors measured at a single site."

dataset_url: "https://www.kaggle.com/competitions/open-problems-multimodal/data"
dataset_reference: lance2024predicting
normalization_methods: [log_cp10k, sqrt_cp10k, l1_sqrt]
output_mod1: '$id/dataset_mod1.h5ad'
output_mod2: '$id/dataset_mod2.h5ad'
output_meta_mod1: '$id/dataset_metadata_mod1.yaml'
output_meta_mod2: '$id/dataset_metadata_mod2.yaml'
output_state: '$id/state.yaml'
publish_dir: s3://openproblems-data/resources/datasets/multimodal
HERE

cat > /tmp/nextflow.config << HERE
process {
  withName:'.*publishStatesProc' {
    memory = '16GB'
    disk = '100GB'
  }
}
HERE

tw launch https://github.com/openproblems-bio/openproblems.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/scrnaseq/op3_loader/main.nf \
  --workspace 53907369739130 \
  --params-file "$params_file" \
  --config /tmp/nextflow.config \
  --labels op3_loader,dataset_loader \
