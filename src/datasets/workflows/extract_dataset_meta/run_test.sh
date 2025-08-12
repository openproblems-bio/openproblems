#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# export TOWER_WORKSPACE_ID=53907369739130

OUTPUT_DIR="output/temp"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

DATASETS_DIR="resources_test/common/pancreas/dataset.h5ad"

export NXF_VER=22.04.5
nextflow run . \
  -main-script target/nextflow/datasets/workflows/extract_dataset_meta/main.nf \
  -profile docker \
  -resume \
  -c src/wf_utils/labels_ci.config \
  --input $DATASETS_DIR \
  --output meta.yaml \
  --publish_dir "$OUTPUT_DIR"