#!/bin/bash

DATASET_DIR="resources_test/common"

#make sure the following command has been executed
#viash ns build -q 'datasets|common' --parallel --setup cb

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# download full dataset as temp file
mkdir -p "$DATASET_DIR/neurips2021_bmmc_cite"

INPUT="$DATASET_DIR/neurips2021_bmmc_cite/temp_neurips2021_bmmc_cite.h5ad"
INPUT_URL="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE194122&format=file&file=GSE194122%5Fopenproblems%5Fneurips2021%5Fcite%5FBMMC%5Fprocessed%2Eh5ad%2Egz"
if [ ! -f "$INPUT" ]; then
  echo "Downloading neurips2021_bmmc_cite dataset"
  
  wget "$INPUT_URL" -O "${INPUT}.gz"
    
  gunzip "${INPUT}.gz"
fi


# download dataset
nextflow run . \
  -main-script target/nextflow/datasets/workflows/process_openproblems_neurips2021_bmmc/main.nf \
  -profile docker \
  -c src/wf_utils/labels_ci.config \
  -resume \
  --id neurips2021_bmmc_cite \
  --input "$INPUT" \
  --mod1 "GEX" \
  --mod2 "ADT" \
  --dataset_name "bmcc (CITE-Seq)" \
  --dataset_url "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE194122" \
  --dataset_reference "Neurips" \
  --dataset_summary "neurips small summary" \
  --dataset_description "neurips big description" \
  --do_subsample true \
  --n_obs 600 \
  --n_vars 1500 \
  --seed 123 \
  --normalization_methods log_cp10k \
  --output_rna '$id/dataset_rna.h5ad' \
  --output_other_mod '$id/dataset_other_mod.h5ad' \
  --output_meta_rna '$id/dataset_metadata_rna.yaml' \
  --output_meta_other_mod '$id/dataset_metadata_other_mod.yaml' \
  --output_state '$id/state.yaml' \
  --publish_dir "$DATASET_DIR"

# run task process dataset components
src/tasks/predict_modality/resources_test_scripts/neurips2021_bmmc.sh