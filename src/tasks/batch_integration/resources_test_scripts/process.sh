#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

RAW_DATA=resources_test/common
DATASET_DIR=resources_test/batch_integration

mkdir -p $DATASET_DIR

# process dataset
echo Running process_dataset
nextflow run . \
  -main-script target/nextflow/batch_integration/workflows/process_datasets/main.nf \
  -profile docker \
  -entry auto \
  --input_states "$RAW_DATA/**/state.yaml" \
  --rename_keys 'input:output_dataset' \
  --settings '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad"}' \
  --publish_dir "$DATASET_DIR" \
  --output_state '$id/state.yaml'
# output_state should be moved to settings once workaround is solved

for id in pancreas cxg_mouse_pancreas_atlas; do
  if [ ! -f $DATASET_DIR/$id/dataset.h5ad ]; then
    echo "Dataset $id not found"
    exit 1
  fi

  echo Running BBKNN on $id
  viash run src/tasks/batch_integration/methods/bbknn/config.vsh.yaml -- \
    --input $DATASET_DIR/$id/dataset.h5ad \
    --output $DATASET_DIR/$id/integrated_graph.h5ad

  echo Running SCVI on $id
  viash run src/tasks/batch_integration/methods/scvi/config.vsh.yaml -- \
    --input $DATASET_DIR/$id/dataset.h5ad \
    --output $DATASET_DIR/$id/integrated_embedding.h5ad

  echo Running combat on $id
  viash run src/tasks/batch_integration/methods/combat/config.vsh.yaml -- \
    --input $DATASET_DIR/$id/dataset.h5ad \
    --output $DATASET_DIR/$id/integrated_feature.h5ad
done