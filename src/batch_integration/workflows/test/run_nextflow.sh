#!/bin/bash

# Run this prior to executing this script:
# bin/viash_build -q 'batch_integration'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# choose a particular version of nextflow
export NXF_VER=21.10.6

bin/nextflow \
  run . \
  -main-script src/batch_integration/workflows/run/main.nf \
  --publishDir output/batch_integration \
  -resume \
  -with-docker

