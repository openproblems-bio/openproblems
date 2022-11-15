#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=22.04.5

bin/nextflow \
  run . \
  -main-script src/datasets/workflows/process_openproblems_v1/main.nf \
  -resume \
  -profile docker \
  -params-file src/datasets/loaders/openproblems_v1/datasets.yaml \
  --publish_dir output/datasets