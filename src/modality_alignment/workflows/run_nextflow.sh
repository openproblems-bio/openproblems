#!/bin/bash

# Run this prior to executing this script:
# bin/project_build -q 'modality_alignment|utils'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

NXF_VER=20.10.0 nextflow run src/modality_alignment/workflows/main.nf \
  -resume \
  --output output/modality_alignment

