#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

DATASET_DIR=resources_test/common/multimodal

set -e

mkdir -p $DATASET_DIR

# download dataset
nextflow run . \
  -main-script target/nextflow/datasets/workflows/process_openproblems_v1_multimodal/main.nf \
  -profile docker \
  -resume \
  --id scicar_cell_lines \
  --input_id scicar_cell_lines \
  --obs_tissue "source" \
  --layer_counts "counts" \
  --obs_cell_type "cell_name" \
  --var_feature_id "index" \
  --var_feature_name "gene_short_name" \
  --dataset_name "sci-CAR cell lines" \
  --dataset_url "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117089" \
  --dataset_reference "cao2018joint" \
  --dataset_summary "sciCAR is a combinatorial indexing-based assay that jointly measures cellular transcriptomes and the accessibility of cellular chromatin in the same cells" \
  --dataset_description "sciCAR is a combinatorial indexing-based assay that jointly measures cellular transcriptomes and the accessibility of cellular chromatin in the same cells. Here, we use two sciCAR datasets that were obtained from the same study. The first dataset contains 4,825 cells from three cell lines (HEK293T cells, NIH/3T3 cells, and A549 cells) at multiple timepoints (0, 1 hour, 3 hours) after dexamethasone treatment. The second dataset contains 11,233 cells from wild-type adult mouse kidney." \
  --dataset_organism "[homo_sapiens, mus_musculus]" \
  --mod1 GEX \
  --mod2 ATAC \
  --do_subsample true \
  --n_obs 600 \
  --n_vars 1500 \
  --seed 123 \
  --normalization_methods log_cp10k \
  --output_mod1 '$id/dataset_mod1.h5ad' \
  --output_mod2 '$id/dataset_mod2.h5ad' \
  --output_meta_mod1 '$id/dataset_metadata_mod1.yaml' \
  --output_meta_mod2 '$id/dataset_metadata_mod2.yaml' \
  --output_state '$id/state.yaml' \
  --publish_dir "$DATASET_DIR"

# run task process dataset components
src/tasks/match_modalities/resources_test_scripts/scicar_cell_lines.sh