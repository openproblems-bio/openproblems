#!/bin/bash

# fail on error
set -e

# ensure we're in the root of the repo
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

DATASET_DIR="s3://openproblems-data/resources/datasets/"

for LOADER in $(aws s3 ls $DATASET_DIR); do

  if [ "$LOADER" == "PRE" ]; then
    continue
  fi

  BASE_DIR="${DATASET_DIR%/}/$LOADER"

  for DATASET in $(aws s3 ls $BASE_DIR); do
    
    if [ "$DATASET" == "PRE" ]; then
      continue
    fi
    
    FILE_DIR="${BASE_DIR%/}/${DATASET%/}/log_cp10k/"
    FILES=$(aws s3 ls $FILE_DIR)
    metafiles=$(echo "$FILES" | grep "meta" | awk '{print $NF}')
    # metafiles=$(find $INPUT -type f -name "*meta*")
    # echo $metafiles

    for metafile in $metafiles; do
      INPUT="${FILE_DIR%/}/$metafile"
      OUTPUT_DIR="../website/datasets/$LOADER/${DATASET%/}/data/"
      OUTPUT_FILE="${metafile%.*}.json"
      echo "Processing $LOADER - $DATASET : $INPUT"

      # start the 
      NXF_VER=23.10.0 nextflow run . \
      -main-script target/nextflow/reporting/process_dataset_metadata/main.nf \
      -profile docker \
      -c src/wf_utils/labels_ci.config \
      --id "extract_metadata" \
      --input "$INPUT" \
      --output "$OUTPUT_FILE" \
      --output_state "state.yaml" \
      --publish_dir "$OUTPUT_DIR"
    done

# cause quarto rerender to index page when in preview mode
# touch ../website/results/$TASK/index.qmd
    done
done