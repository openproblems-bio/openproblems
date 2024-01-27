#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

DATASET_DIR=resources_test/common

set -e

mkdir -p $DATASET_DIR

wget https://raw.githubusercontent.com/theislab/scib/c993ffd9ccc84ae0b1681928722ed21985fb91d1/scib/resources/g2m_genes_tirosh_hm.txt -O $DATASET_DIR/temp_g2m_genes_tirosh_hm.txt
wget https://raw.githubusercontent.com/theislab/scib/c993ffd9ccc84ae0b1681928722ed21985fb91d1/scib/resources/s_genes_tirosh_hm.txt -O $DATASET_DIR/temp_s_genes_tirosh_hm.txt
KEEP_FEATURES=`cat $DATASET_DIR/temp_g2m_genes_tirosh_hm.txt $DATASET_DIR/temp_s_genes_tirosh_hm.txt | paste -sd ":" -`

# download dataset
nextflow run . \
  -main-script target/nextflow/datasets/workflows/process_openproblems_v1/main.nf \
  -profile docker \
  -resume \
  --id pancreas \
  --input_id pancreas \
  --obs_cell_type "celltype" \
  --obs_batch "tech" \
  --layer_counts "counts" \
  --dataset_name "Human pancreas" \
  --dataset_url "https://theislab.github.io/scib-reproducibility/dataset_pancreas.html" \
  --dataset_reference "luecken2022benchmarking" \
  --dataset_summary "Human pancreas cells dataset from the scIB benchmarks" \
  --dataset_description "Human pancreatic islet scRNA-seq data from 6 datasets across technologies (CEL-seq, CEL-seq2, Smart-seq2, inDrop, Fluidigm C1, and SMARTER-seq)." \
  --dataset_organism "homo_sapiens" \
  --keep_cell_type_categories "acinar:beta" \
  --keep_batch_categories "celseq:inDrop4:smarter" \
  --keep_features "$KEEP_FEATURES" \
  --seed 123 \
  --normalization_methods log_cp10k \
  --do_subsample true \
  --output_raw '$id/raw.h5ad' \
  --output_normalized '$id/normalized.h5ad' \
  --output_hvg '$id/hvg.h5ad' \
  --output_pca '$id/pca.h5ad' \
  --output_knn '$id/knn.h5ad' \
  --output_dataset '$id/dataset.h5ad' \
  --output_meta '$id/dataset_meta.yaml' \
  --output_state '$id/state.yaml' \
  --publish_dir "$DATASET_DIR"

rm -r $DATASET_DIR/temp_*

# run task process dataset components
src/tasks/batch_integration/resources_test_scripts/process.sh
src/tasks/denoising/resources_test_scripts/pancreas.sh
src/tasks/dimensionality_reduction/resources_test_scripts/pancreas.sh
src/tasks/label_projection/resources_test_scripts/pancreas.sh
