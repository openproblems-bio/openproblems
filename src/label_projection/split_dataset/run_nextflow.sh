#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=22.04.5

bin/nextflow \
  run . \
  -main-script target/nextflow/label_projection/split_dataset/main.nf \
  -profile docker \
  -resume \
  -params-file src/label_projection/split_dataset/params.yaml \
  --publish_dir resources/label_projection/openproblems_v1