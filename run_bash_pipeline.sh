#!/bin/bash

# Run 'bin/project_build' prior to executing this script

TARGET=target/docker/modality_alignment
OUTPUT=out_bash/modality_alignment

mkdir -p $OUTPUT/datasets
mkdir -p $OUTPUT/methods
mkdir -p $OUTPUT/metrics

if [ ! -f "$OUTPUT/dataset/citeseq_cbmc.h5ad" ]; then
  "$TARGET/datasets/citeseq_cbmc/citeseq_cbmc" --output "$OUTPUT/dataset/citeseq_cbmc.h5ad"
fi


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


for met in `ls $TARGET/metrics`; do
  for outp in `ls $OUTPUT/methods`; do
    out_id="${outp%.*}"
    input_h5ad="$OUTPUT/methods/$out_id.h5ad"
    output_h5ad="$OUTPUT/metrics/${out_id}_$met.h5ad"
    if [ ! -f "$output_h5ad" ]; then
      echo "> $TARGET/metric/$met/$met" -i "$input_h5ad" -o "$output_h5ad"
      "$TARGET/metric/$met/$met" -i "$input_h5ad" -o "$output_h5ad"
    fi
  done
done
