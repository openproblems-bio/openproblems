#!/bin/bash
#It is an example how to set parameters and execute processing for the op3 dataset.
set -e

params_file="/tmp/datasets_op3.yaml"

cat > "$params_file" << 'HERE'
param_list:
  - id: op3
    input: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279945/suppl/GSE279945_sc_counts_processed.h5ad
    dataset_name: "OP3: single-cell multimodal dataset in PBMCs for perturbation prediction benchmarking"
    dataset_summary: "The Open Problems Perurbation Prediction (OP3) dataset with small molecule perturbations in PBMCs"
    dataset_description: "The OP3 dataset is to-date the largest single-cell small molecule perturbation dataset in primary tissue with multiple donor replicates."
    dataset_url: "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279945/suppl/GSE279945_sc_counts_processed.h5ad"
    dataset_reference: GSE279945
normalization_methods: 
  - log_cp10k
output_dataset: '$id/dataset.h5ad'
do_subsample: False
output_meta: '$id/dataset_metadata.yaml'
output_state: '$id/state.yaml'
output_raw: force_null
output_normalized: force_null
output_pca: force_null
output_hvg: force_null
output_knn: force_null
publish_dir: s3://openproblems-data/resources/datasets/op3
HERE

cat > "/tmp/nextflow.config" << 'HERE'
process {
  withName:'.*publishStatesProc' {
    memory = '100GB'
    disk = '1000GB'
  }
}
HERE

tw launch https://github.com/openproblems-bio/openproblems.git \
 --revision main_build \
 --pull-latest \
 --main-script target/nextflow/datasets/workflows/scrnaseq/process_op3/main.nf \
 --workspace 53907369739130 \
 --params-file "$params_file" \
 --labels op3,dataset_loader \
 --config /tmp/nextflow.config

# set -x
# nextflow run . \
#   -main-script target/nextflow/datasets/workflows/scrnaseq/process_op3/main.nf \
#   -profile docker \
#   -resume \
#   -params-file "$params_file" \
#   -config /tmp/nextflow.config
