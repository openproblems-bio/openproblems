#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export TOWER_WORKSPACE_ID=53907369739130

OUTPUT_DIR="resources/datasets/multimodal"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

params_file="/tmp/datasets_openproblems_v1_multimodal_params.yaml"

cat > "$params_file" << 'HERE'
param_list:
  - id: openproblems_v1_multimodal/citeseq_cbmc
    dataset_name: "CITE-Seq CBMC"
    dataset_summary: "CITE-seq profiles of 8k Cord Blood Mononuclear Cells"
    dataset_description: "8k cord blood mononuclear cells profiled by CITEseq using a panel of 13 antibodies."
    dataset_reference: stoeckius2017simultaneous
    dataset_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100866
    dataset_organism: homo_sapiens
    layer_counts: counts

normalization_methods: [log_cp10k, sqrt_cp10k, l1_sqrt]
output_mod1: '$id/dataset_mod1.h5ad'
output_mod2: '$id/dataset_mod2.h5ad'
output_meta_mod1: '$id/dataset_metadata_mod1.yaml'
output_meta_mod2: '$id/dataset_metadata_mod2.yaml'
output_state: '$id/state.yaml'
HERE

export NXF_VER=22.04.5
nextflow \
  run . \
  -main-script target/nextflow/datasets/workflows/multimodal/process_openproblems_v1_multimodal/main.nf \
  -profile docker \
  -resume \
  -params-file "$params_file" \
  --publish_dir "$OUTPUT_DIR"
