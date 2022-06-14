#!/bin/bash
#
#make sure the following command has been executed
#bin/viash_build -q 'label_projection|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

target/docker/common/data_loader/data_loader\
    --url "https://ndownloader.figshare.com/files/24539828"\
    --name "pancreas"\
    --output src/label_projection/resources/raw_data.h5ad
