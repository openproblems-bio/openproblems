#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export TOWER_WORKSPACE_ID=53907369739130

OUTPUT_DIR="resources/datasets/openproblems_v1"

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

params_file="$OUTPUT_DIR/params.yaml"

if [ ! -f $params_file ]; then
  cat > "$params_file" << 'HERE'
param_list:
  - id: allen_brain_atlas
    obs_celltype: label
    layer_counts: counts
    dataset_id: allen_brain_atlas
    dataset_name: Mouse Brain Atlas
    data_url: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71585
    data_reference: tasic2016adult
    dataset_summary: Adult mouse primary visual cortex
    dataset_description: A murine brain atlas with adjacent cell types as assumed benchmark truth, inferred from deconvolution proportion correlations using matching 10x Visium slides (see Dimitrov et al., 2022).
    dataset_organism: mus_musculus

  - id: cengen
    obs_celltype: cell_type
    obs_batch: experiment_code
    obs_tissue: tissue
    layer_counts: counts
    dataset_id: cengen
    dataset_name: CeNGEN
    data_url: https://www.cengen.org
    data_reference: hammarlund2018cengen
    dataset_summary: Complete Gene Expression Map of an Entire Nervous System
    dataset_description: 100k FACS-isolated C. elegans neurons from 17 experiments sequenced on 10x Genomics.
    dataset_organism: caenorhabditis_elegans

  - id: immune_cells
    obs_celltype: final_annotation
    obs_batch: batch
    obs_tissue: tissue
    layer_counts: counts
    dataset_id: immune_cells
    dataset_name: Human immune
    data_url: https://theislab.github.io/scib-reproducibility/dataset_immune_cell_hum.html
    data_reference: luecken2022benchmarking
    dataset_summary: Human immune cells dataset from the scIB benchmarks
    dataset_description: Human immune cells from peripheral blood and bone marrow taken from 5 datasets comprising 10 batches across technologies (10X, Smart-seq2).
    dataset_organism: homo_sapiens

  - id: mouse_blood_olsson_labelled
    obs_celltype: celltype
    layer_counts: counts
    dataset_id: mouse_blood_olsson_labelled
    dataset_name: Mouse myeloid
    data_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70245
    data_reference: olsson2016single
    dataset_summary: Myeloid lineage differentiation from mouse blood
    dataset_description: 660 FACS-isolated myeloid cells from 9 experiments sequenced using C1 Fluidigm and SMARTseq in 2016 by Olsson et al.
    dataset_organism: mus_musculus

  - id: mouse_hspc_nestorowa2016
    obs_celltype: cell_type_label
    layer_counts: counts
    dataset_id: mouse_hspc_nestorowa2016
    dataset_name: Mouse HSPC
    data_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81682
    data_reference: nestorowa2016single
    dataset_summary: Haematopoeitic stem and progenitor cells from mouse bone marrow
    dataset_description: 1656 hematopoietic stem and progenitor cells from mouse bone marrow. Sequenced by Smart-seq2. 
    dataset_organism: mus_musculus

  - id: pancreas
    obs_celltype: celltype
    obs_batch: tech
    layer_counts: counts
    dataset_id: pancreas
    dataset_name: Human pancreas
    data_url: https://theislab.github.io/scib-reproducibility/dataset_pancreas.html
    data_reference: luecken2022benchmarking
    dataset_summary: Human pancreas cells dataset from the scIB benchmarks
    dataset_description: Human pancreatic islet scRNA-seq data from 6 datasets across technologies (CEL-seq, CEL-seq2, Smart-seq2, inDrop, Fluidigm C1, and SMARTER-seq). 
    dataset_organism: homo_sapiens

  - id: tabula_muris_senis_droplet_lung
    obs_celltype: cell_type
    obs_batch: donor_id
    layer_counts: counts
    dataset_id: tabula_muris_senis_droplet_lung
    dataset_name: Tabula Muris Senis Lung
    data_url: https://tabula-muris-senis.ds.czbiohub.org
    data_reference: tabula2020single
    dataset_summary: Aging mouse lung cells from Tabula Muris Senis
    dataset_description: All lung cells from 10x profiles in Tabula Muris Senis, a 500k cell-atlas from 18 organs and tissues across the mouse lifespan.
    dataset_organism: mus_musculus

  - id: tenx_1k_pbmc
    layer_counts: counts
    dataset_id: tenx_1k_pbmc
    dataset_name: 1k PBMCs
    data_url: https://www.10xgenomics.com/resources/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0
    data_reference: 10x2018pbmc
    dataset_summary: 1k peripheral blood mononuclear cells from a healthy donor
    dataset_description: 1k Peripheral Blood Mononuclear Cells (PBMCs) from a healthy donor. Sequenced on 10X v3 chemistry in November 2018 by 10X Genomics.
    dataset_organism: homo_sapiens

  - id: tenx_5k_pbmc
    layer_counts: counts
    dataset_id: tenx_5k_pbmc
    dataset_name: 5k PBMCs
    data_url: https://www.10xgenomics.com/resources/datasets/5-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-with-cell-surface-proteins-v-3-chemistry-3-1-standard-3-1-0
    data_reference: 10x2019pbmc
    dataset_summary: 5k peripheral blood mononuclear cells from a healthy donor
    dataset_description: 5k Peripheral Blood Mononuclear Cells (PBMCs) from a healthy donor. Sequenced on 10X v3 chemistry in July 2019 by 10X Genomics.
    dataset_organism: homo_sapiens

  - id: tnbc_wu2021
    obs_celltype: celltype_minor
    layer_counts: counts
    dataset_id: tnbc_wu2021
    dataset_name: Triple-Negative Breast Cancer
    data_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118389
    data_reference: wu2021single
    dataset_summary: 1535 cells from six fresh triple-negative breast cancer tumors.
    dataset_description: 1535 cells from six TNBC donors by (Wu et al., 2021). This dataset includes cytokine activities, inferred using a multivariate linear model with cytokine-focused signatures, as assumed true cell-cell communication (Dimitrov et al., 2022).
    dataset_organism: homo_sapiens
    
  - id: zebrafish
    obs_celltype: cell_type
    obs_batch: lab
    layer_counts: counts
    dataset_id: zebrafish
    dataset_name: Zebrafish embryonic cells
    data_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112294
    data_reference: wagner2018single
    dataset_summary: Single-cell mRNA sequencing of zebrafish embryonic cells.
    dataset_description: 90k cells from zebrafish embryos throughout the first day of development, with and without a knockout of chordin, an important developmental gene. 
    dataset_organism: danio_rerio

output: '$id.h5ad'
HERE
fi

export NXF_VER=22.04.5
nextflow \
  run . \
  -main-script src/datasets/workflows/process_openproblems_v1/main.nf \
  -profile docker \
  -resume \
  -params-file "$params_file" \
  --publish_dir "$OUTPUT_DIR"
  
  # -with-tower
