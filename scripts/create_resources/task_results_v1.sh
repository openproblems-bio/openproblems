#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

OUT_DIR="resources_test/openproblems/task_results_v1"

set -e

echo ">> Removing previously generated results"
if [ -d "$OUT_DIR" ]; then
  rm -r "$OUT_DIR"
fi

echo ">> Fetch results in v1 format"
mkdir -p "$OUT_DIR/data/"
TMPDIR=$(mktemp -d)

wget https://github.com/openproblems-bio/website/archive/refs/tags/v2.3.6.zip -O "$TMPDIR/website-v2.3.6.zip"
unzip "$TMPDIR/website-v2.3.6.zip" -d "$TMPDIR"
cp -r "$TMPDIR/website-2.3.6/results/batch_integration_embed/data/" "$OUT_DIR/processed"

echo ">> Uploading results to S3"
aws s3 sync --profile op \
  "resources_test/openproblems/task_results_v1/" \
  "s3://openproblems-data/resources_test/openproblems/task_results_v1/" \
  --delete --dryrun
