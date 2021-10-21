#!/bin/bash

# Run this prior to executing this script:
# bin/viash_build -q 'modality_alignment|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# choose a particular version of nextflow
export NXF_VER=21.04.1

bin/nextflow \
  run https://github.com/openproblems-bio/opsca-viash.git \
  -r main_build \
  -main-script src/modality_alignment/workflows/main.nf \
  -c src/modality_alignment/workflows/nextflow.config \
  --output output/modality_alignment \
  -resume

