#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export TOWER_WORKSPACE_ID=53907369739130

OUTPUT_DIR="resources/datasets/openproblems_v1"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

params_file="$OUTPUT_DIR/params.yaml"

if [ ! -f $params_file ]; then
  cat > "$params_file" << 'HERE'
param_list:
  - id: allen_brain_atlas
    obs_celltype: label
    layer_counts: counts

  - id: cengen
    obs_celltype: cell_type
    obs_batch: experiment_code
    obs_tissue: tissue
    layer_counts: counts

  - id: immune_cells
    obs_celltype: final_annotation
    obs_batch: batch
    obs_tissue: tissue
    layer_counts: counts

  - id: mouse_blood_olssen_labelled
    obs_celltype: celltype
    layer_counts: counts

  - id: mouse_hspc_nestorowa2016
    obs_celltype: cell_type_label
    layer_counts: counts

  - id: pancreas
    obs_celltype: celltype
    obs_batch: tech
    layer_counts: counts

  - id: tabula_muris_senis_droplet_lung
    obs_celltype: cell_type
    obs_batch: donor_id
    layer_counts: counts

  - id: tenx_1k_pbmc
    layer_counts: counts

  - id: tenx_5k_pbmc
    layer_counts: counts

  - id: tnbc_wu2021
    obs_celltype: celltype_minor
    layer_counts: counts
    
  - id: zebrafish
    obs_celltype: cell_type
    obs_batch: lab
    layer_counts: counts

output: '$id.h5ad'
HERE
fi

export NXF_VER=22.04.5
nextflow \
  run . \
  -main-script src/datasets/workflows/process_openproblems_v1/main.nf \
  -profile docker \
  -resume \
  -params-file "$params_file" \
  --publish_dir "$OUTPUT_DIR" \
  -with-tower
