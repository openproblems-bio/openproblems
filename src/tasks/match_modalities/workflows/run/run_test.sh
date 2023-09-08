#!/bin/bash

# Run this prior to executing this script:
# viash ns build -q 'match_modalities|common' --setup cb --parallel

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

DATASET_DIR=resources_test/common/multimodal

# choose a particular version of nextflow
export NXF_VER=23.04.2

nextflow \
  run . \
  -main-script src/tasks/match_modalities/workflows/run/main.nf \
  -resume \
  -c src/wf_utils/labels_ci.config \
  -profile docker \
  --id scicar \
  --input_mod1 $DATASET_DIR/dataset_mod1.h5ad \
  --input_mod2 $DATASET_DIR/dataset_mod2.h5ad \
  --output scores.tsv \
  --publish_dir output/match_modalities/ \



