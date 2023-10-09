#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export TOWER_WORKSPACE_ID=53907369739130

OUTPUT_DIR="resources/datasets/openproblems_v1_multimodal"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

params_file="$OUTPUT_DIR/params.yaml"

if [ ! -f $params_file ]; then
  cat > "$params_file" << 'HERE'
param_list:
  - id: citeseq_cbmc
    dataset_id: citeseq_cbmc
    dataset_name: "CITE-Seq CBMC"
    dataset_summary: "CITE-seq profiles of 8k Cord Blood Mononuclear Cells"
    dataset_description: "8k cord blood mononuclear cells profiled by CITEsequsing a panel of 13 antibodies."
    data_reference: stoeckius2017simultaneous
    data_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100866
    dataset_organism: homo_sapiens
    layer_counts: counts

  - id: scicar_cell_lines
    dataset_id: scicar_cell_lines
    dataset_name: "sci-CAR Cell Lines"
    dataset_summary: "sci-CAR profiles of 5k cell line cells (HEK293T, NIH/3T3, A549) across three treatment conditions (DEX 0h, 1h and 3h)"
    dataset_description: "Single cell RNA-seq and ATAC-seq co-profiling for HEK293T cells, NIH/3T3 cells, A549 cells across three treatment conditions (DEX 0 hour, 1 hour and 3 hour treatment)."
    data_reference: cao2018joint
    data_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117089
    dataset_organism: [homo_sapiens, mus_musculus]
    obs_celltype: cell_name
    layer_counts: counts

  - id: scicar_mouse_kidney
    dataset_id: scicar_mouse_kidney
    dataset_name: "sci-CAR Mouse Kidney"
    dataset_summary: "sci-CAR profiles of 11k mouse kidney cells"
    dataset_description: "Single cell RNA-seq and ATAC-seq co-profiling of 11k mouse kidney cells."
    data_reference: cao2018joint
    data_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117089
    dataset_organism: mus_musculus
    obs_celltype: cell_name
    obs_batch: replicate
    layer_counts: counts

normalization_id: [log_cp10k, sqrt_cp10k, l1_sqrt]
output_dataset_mod1: '$id/dataset_mod1.h5ad'
output_dataset_mod1: '$id/dataset_mod2.h5ad'
output_meta_mod1: '$id/dataset_metadata_mod1.h5ad'
output_meta_mod1: '$id/dataset_metadata_mod2.h5ad'
output_state: '$id/state.yaml'
HERE
fi

export NXF_VER=22.04.5
nextflow \
  run . \
  -main-script target/nextflow/datasets/workflows/process_openproblems_v1_multimodal/main.nf \
  -profile docker \
  -resume \
  -params-file "$params_file" \
  --publish_dir "$OUTPUT_DIR"
