#!/bin/bash

# Run this prior to executing this script:
# bin/project_build -q 'modality_alignment|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

NXF_VER=20.10.0 bin/nextflow run src/trajectory_inference/workflows/main.nf \
  -resume \
  --publish_dir output/trajectory_inference

