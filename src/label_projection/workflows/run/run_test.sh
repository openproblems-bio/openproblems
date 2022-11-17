#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=22.04.5

bin/nextflow \
  run . \
  -main-script src/label_projection/workflows/run/main.nf \
  -profile docker \
  -resume \
  -params-file src/label_projection/workflows/run/params_test.yaml \
  --publish_dir output/label_projection