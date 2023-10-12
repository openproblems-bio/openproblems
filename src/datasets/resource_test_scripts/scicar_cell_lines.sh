#!/bin/bash
#
#make sure the following command has been executed
#viash ns build -q 'datasets|common' --parallel --setup cb

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

DATASET_DIR=resources_test/common

set -e

mkdir -p $DATASET_DIR

# download dataset
nextflow run . \
  -main-script target/nextflow/datasets/workflows/process_openproblems_v1_multimodal/main.nf \
  -profile docker \
  -resume \
  --id scicar_cell_lines \
  --obs_tissue "source" \
  --layer_counts "counts" \
  --obs_celltype "cell_name" \
  --dataset_id scicar_cell_lines \
  --dataset_name "sci-CAR cell lines" \
  --data_url "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117089" \
  --data_reference "cao2018joint" \
  --dataset_summary "sciCAR is a combinatorial indexing-based assay that jointly measures cellular transcriptomes and the accessibility of cellular chromatin in the same cells" \
  --dataset_description "sciCAR is a combinatorial indexing-based assay that jointly measures cellular transcriptomes and the accessibility of cellular chromatin in the same cells. Here, we use two sciCAR datasets that were obtained from the same study. The first dataset contains 4,825 cells from three cell lines (HEK293T cells, NIH/3T3 cells, and A549 cells) at multiple timepoints (0, 1 hour, 3 hours) after dexamethasone treatment. The second dataset contains 11,233 cells from wild-type adult mouse kidney." \
  --dataset_organism "[homo_sapiens, mus_musculus]" \
  --do_subsample true \
  --n_obs 600 \
  --n_vars 1500 \
  --seed 123 \
  --normalization_methods log_cp10k \
  --output_dataset_mod1 '$id/dataset_mod1.h5ad' \
  --output_dataset_mod2 '$id/dataset_mod2.h5ad' \
  --output_meta_mod1 '$id/dataset_metadata_mod1.yaml' \
  --output_meta_mod2 '$id/dataset_metadata_mod2.yaml' \
  --output_state '$id/state.yaml' \
  --publish_dir "$DATASET_DIR"

# run task process dataset components
src/tasks/match_modalities/resources_test_scripts/scicar_cell_lines.sh