#!/bin/bash

# Run this prior to executing this script:
# bin/viash_build -q 'modality_alignment|utils'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# choose a particular version of nextflow
export NXF_VER=21.10.6

bin/nextflow \
  run . \
  -main-script src/label_projection/workflows/run/main.nf \
  --publishDir output/label_projection \
  -resume \
  -with-docker
