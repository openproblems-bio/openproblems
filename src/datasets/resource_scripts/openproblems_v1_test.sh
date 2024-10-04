#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export TOWER_WORKSPACE_ID=53907369739130

OUTPUT_DIR="resources/datasets/scrnasrq"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

params_file="/tmp/datasets_openproblems_v1_params.yaml"

cat > "$params_file" << 'HERE'
param_list:
  - id: openproblems_v1/pancreas
    obs_cell_type: celltype
    obs_batch: tech
    layer_counts: counts
    dataset_name: Human pancreas
    dataset_url: https://theislab.github.io/scib-reproducibility/dataset_pancreas.html
    dataset_reference: luecken2022benchmarking
    dataset_summary: Human pancreas cells dataset from the scIB benchmarks
    dataset_description: Human pancreatic islet scRNA-seq data from 6 datasets across technologies (CEL-seq, CEL-seq2, Smart-seq2, inDrop, Fluidigm C1, and SMARTER-seq). 
    dataset_organism: homo_sapiens

normalization_methods: [log_cp10k, sqrt_cp10k, l1_sqrt]
output_dataset: '$id/dataset.h5ad'
output_meta: '$id/dataset_metadata.yaml'
output_state: '$id/state.yaml'
output_raw: force_null
output_normalized: force_null
output_pca: force_null
output_hvg: force_null
output_knn: force_null
HERE

export NXF_VER=23.04.2
nextflow run . \
  -main-script target/nextflow/datasets/workflows/process_openproblems_v1/main.nf \
  -profile docker \
  -resume \
  -params-file "$params_file" \
  --publish_dir "$OUTPUT_DIR"
  
  # -with-tower
