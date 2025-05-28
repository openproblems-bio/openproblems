#!/bin/bash

REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

DATASET_DIR=op3_data/

set -e
mkdir -p $DATASET_DIR

set -x
nextflow run . \
  -main-script target/nextflow/datasets/workflows/scrnaseq/process_op3/main.nf \
  -profile docker \
  -resume \
  --input https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279945/suppl/GSE279945_sc_counts_processed.h5ad \
  --id op3 \
  --dataset_name "OP3: single-cell multimodal dataset in PBMCs for perturbation prediction benchmarking" \
  --dataset_summary: "The Open Problems Perurbation Prediction (OP3) dataset with small molecule perturbations in PBMCs" \
  --dataset_description: "The OP3 dataset is to-date the largest single-cell small molecule perturbation dataset in primary tissue with multiple donor replicates." \
  --dataset_url: "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279945/suppl/GSE279945_sc_counts_processed.h5ad" \
  --dataset_reference: GSE279945 \
  --normalization_methods: log_cp10k \
  --output_dataset: '$id/dataset.h5ad' \
  --output_meta: '$id/dataset_metadata.yaml' \
  --output_state: '$id/state.yaml' \
  --output_raw: force_null \
  --output_normalized: force_null \
  --output_pca: force_null \
  --output_hvg: force_null \
  --output_knn: force_null \
  --publish_dir: '$DATASET_DIR' \
  -config common/nextflow_helpers/labels_ci.config

