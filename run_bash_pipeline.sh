#!/bin/bash

# Run this prior to executing this script:
# bin/project_build

TARGET=target/docker/modality_alignment
OUTPUT=out_bash/modality_alignment

mkdir -p $OUTPUT/datasets
mkdir -p $OUTPUT/methods
mkdir -p $OUTPUT/metrics

# generate datasets
if [ ! -f "$OUTPUT/datasets/citeseq_cbmc.h5ad" ]; then
  "$TARGET/datasets/citeseq_cbmc/citeseq_cbmc" --output "$OUTPUT/datasets/citeseq_cbmc.h5ad"
fi

# run all methods on all datasets
for meth in `ls "$TARGET/methods"`; do
  for dat in `ls "$OUTPUT/datasets"`; do
    dat_id="${dat%.*}"
    input_h5ad="$OUTPUT/datasets/$dat_id.h5ad"
    output_h5ad="$OUTPUT/methods/${dat_id}_$meth.h5ad"
    if [ ! -f "$output_h5ad" ]; then
      echo "> $TARGET/methods/$meth/$meth -i $input_h5ad -o $output_h5ad"
      "$TARGET/methods/$meth/$meth" -i "$input_h5ad" -o "$output_h5ad"
    fi
  done
done

# run all metrics on all outputs
for met in `ls "$TARGET/metrics"`; do
  for outp in `ls "$OUTPUT/methods"`; do
    out_id="${outp%.*}"
    input_h5ad="$OUTPUT/methods/$out_id.h5ad"
    output_h5ad="$OUTPUT/metrics/${out_id}_$met.h5ad"
    if [ ! -f "$output_h5ad" ]; then
      echo "> $TARGET/metrics/$met/$met" -i "$input_h5ad" -o "$output_h5ad"
      "$TARGET/metrics/$met/$met" -i "$input_h5ad" -o "$output_h5ad"
    fi
  done
done

# concatenate all scores into one tsv
INPUTS=$(ls -1 "$OUTPUT/metrics" | sed "s#.*#-i '$OUTPUT/metrics/&'#" | tr '\n' ' ')
eval "$TARGET/utils/docker/utils/extract_scores" $INPUTS -o "$OUTPUT/scores.tsv"