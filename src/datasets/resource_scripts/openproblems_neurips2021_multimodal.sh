#!/bin/bash

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE194122&format=file&file=GSE194122%5Fopenproblems%5Fneurips2021%5Fcite%5FBMMC%5Fprocessed%2Eh5ad%2Egz" \
  -O "/tmp/neurips2021_bmmc_cite.h5ad.gz"

gunzip "/tmp/neurips2021_bmmc_cite.h5ad.gz"

wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE194122&format=file&file=GSE194122%5Fopenproblems%5Fneurips2021%5Fmultiome%5FBMMC%5Fprocessed%2Eh5ad%2Egz" \
  -O "/tmp/neurips2021_bmmc_multiome.h5ad.gz"

gunzip "/tmp/neurips2021_bmmc_multiome.h5ad.gz"



params_file="/tmp/datasets_openproblems_nuerips2021_params.yaml"



cat > "$params_file" << HERE
param_list:
  - id: openproblems_neurips2021/bmmc_cite
    input: "/tmp/neurips2021_bmmc_cite.h5ad"
    mod1: GEX
    mod2: ADT
    dataset_name: bmmc (CITE-Seq)
    dataset_organism: homo_sapiens
    dataset_summary: "Short Summary."
    dataset_description: "Full description."

  - id: openproblems_neurips2021/bmmc_multiome
    input: "/tmp/neurips2021_bmmc_multiome.h5ad"
    mod1: GEX
    mod2: ATAC
    dataset_name: bmmc (Multiome)
    dataset_organism: homo_sapiens
    dataset_summary: "Short Summary."
    dataset_description: "Full description."

dataset_url: "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122"
dataset_reference: neurips
output_rna: '$id/dataset_rna.had'
output_other_mod: '$id/dataset_other_mod.h5ad'
output_meta_rna: '$id/dataset_metadata_rna.yaml'
output_meta_other_mod: '$id/dataset_metadata_other_mod.yaml'
output_state: '$id/state.yaml
HERE

export NXF_VER=23.04.2
nextflow run . \
  -main-script target/nextflow/datasets/workflows/process_openproblems_neurips2021_bmmc/main.nf \
  -profile docker \
  -resume \
  -params-file "$params_file" \
  --publish_dir "resources/datasets/openproblems_neurips2021"

