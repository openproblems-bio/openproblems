#!/bin/bash
#make sure the following command has been executed
#bin/viash_build -q 'label_projection|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

sh src/label_projection/workflows/test/toy_test_data.sh

target/docker/label_projection/data/preprocess/preprocess\
    --input src/label_projection/resources/toy_data.h5ad\
    --output src/label_projection/resources/toy_preprocessed_data.h5ad
