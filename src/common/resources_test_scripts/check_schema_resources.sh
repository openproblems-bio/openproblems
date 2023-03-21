#!/bin/bash

# make sure folloewing command has been executed
# viash ns build -q 'common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

OUTPUT_DIR="resources_test/common/check_schema"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

# Create small git sha input file
correct_file="$OUTPUT_DIR/anndata_correct.yaml"
error_file="$OUTPUT_DIR/anndata_error.yaml"

cat <<EOT > $correct_file
type: file
description: "A preprocessed dataset"
example: "preprocessed.h5ad"
info:
  short_description: "Preprocessed dataset"
  slots:
    layers: 
      - type: integer
        name: counts
        description: Raw counts
    uns:
      - type: string
        name: dataset_id
        description: "A unique identifier for the dataset"
EOT

cat <<EOT > $error_file
type: file
description: "A preprocessed dataset"
example: "preprocessed.h5ad"
info:
  short_description: "Preprocessed dataset"
  slots:
    layers: 
      - type: integer
        name: counts
        description: Raw counts
    uns:
      - type: string
        name: dataset_id
        description: "A unique identifier for the dataset"
      - type: string
        name: error_test
        description: "A made up uns variable to test if error is picked up"
EOT