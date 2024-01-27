#!/bin/bash

# fail on error
set -e

# ensure we're in the root of the repo
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# TODO: Add multimodal datasets
for LOADER in "cellxgene_census" "openproblems_v1"; do 
  BASE_DIR="s3://openproblems-data/resources/datasets/$LOADER/"

  for DATASET in $(aws s3 ls $BASE_DIR); do

    if [ "$DATASET" == "PRE" ]; then
      continue
    fi
  
    INPUT="${BASE_DIR%/}/${DATASET%/}/log_cp10k/dataset_metadata.yaml"
    OUTPUT_DIR="../website/datasets/$LOADER/$DATASET"

    echo "Processing $LOADER - $DATASET : $INPUT"
  # # temp sync
  # aws s3 sync $INPUT_DIR output/temp

    # start the run
    NXF_VER=23.10.0 nextflow run . \
      -main-script target/nextflow/common/process_dataset_metadata/run/main.nf \
      -profile docker \
      -c src/wf_utils/labels_ci.config \
      --id "process" \
      --input "$INPUT" \
      --output_state "state.yaml" \
      --publish_dir "$OUTPUT_DIR"

# cause quarto rerender to index page when in preview mode
# touch ../website/results/$TASK/index.qmd
    done
done