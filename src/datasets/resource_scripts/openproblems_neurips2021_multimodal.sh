#!/bin/bash

params_file="/tmp/datasets_openproblems_neurips2021_params.yaml"

cat > "$params_file" << 'HERE'
param_list:
  - id: openproblems_neurips2021/bmmc_cite
    # input: "/tmp/neurips2021_bmmc_cite.h5ad"
    input: "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194122/suppl/GSE194122%5Fopenproblems%5Fneurips2021%5Fcite%5FBMMC%5Fprocessed%2Eh5ad%2Egz"
    mod1: GEX
    mod2: ADT
    dataset_name: OpenProblems NeurIPS2021 CITE-Seq
    dataset_organism: homo_sapiens
    dataset_summary: Single-cell CITE-Seq (GEX+ADT) data collected from bone marrow mononuclear cells of 12 healthy human donors.
    dataset_description: "Single-cell CITE-Seq data collected from bone marrow mononuclear cells of 12 healthy human donors using the 10X 3 prime Single-Cell Gene Expression kit with Feature Barcoding in combination with the BioLegend TotalSeq B Universal Human Panel v1.0. The dataset was generated to support Multimodal Single-Cell Data Integration Challenge at NeurIPS 2021. Samples were prepared using a standard protocol at four sites. The resulting data was then annotated to identify cell types and remove doublets. The dataset was designed with a nested batch layout such that some donor samples were measured at multiple sites with some donors measured at a single site."

  - id: openproblems_neurips2021/bmmc_multiome
    # input: "/tmp/neurips2021_bmmc_multiome.h5ad"
    input: "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194122/suppl/GSE194122%5Fopenproblems%5Fneurips2021%5Fmultiome%5FBMMC%5Fprocessed%2Eh5ad%2Egz"
    mod1: GEX
    mod2: ATAC
    dataset_name: OpenProblems NeurIPS2021 Multiome
    dataset_organism: homo_sapiens
    dataset_summary: Single-cell Multiome (GEX+ATAC) data collected from bone marrow mononuclear cells of 12 healthy human donors.
    dataset_description: "Single-cell CITE-Seq data collected from bone marrow mononuclear cells of 12 healthy human donors using the 10X Multiome Gene Expression and Chromatin Accessibility kit. The dataset was generated to support Multimodal Single-Cell Data Integration Challenge at NeurIPS 2021. Samples were prepared using a standard protocol at four sites. The resulting data was then annotated to identify cell types and remove doublets. The dataset was designed with a nested batch layout such that some donor samples were measured at multiple sites with some donors measured at a single site."

dataset_url: "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122"
dataset_reference: luecken2021neurips
normalization_methods: [log_cp10k, sqrt_cp10k, l1_sqrt]
output_rna: '$id/dataset_rna.h5ad'
output_other_mod: '$id/dataset_other_mod.h5ad'
output_meta_rna: '$id/dataset_metadata_rna.yaml'
output_meta_other_mod: '$id/dataset_metadata_other_mod.yaml'
output_state: '$id/state.yaml'
publish_dir: s3://openproblems-data/resources/datasets
HERE

cat > /tmp/nextflow.config << HERE
process {
  withName:'.*publishStatesProc' {
    memory = '16GB'
    disk = '100GB'
  }
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_openproblems_neurips2021_bmmc/main.nf \
  --workspace 53907369739130 \
  --compute-env 1pK56PjjzeraOOC2LDZvN2 \
  --params-file "$params_file" \
  --config /tmp/nextflow.config
