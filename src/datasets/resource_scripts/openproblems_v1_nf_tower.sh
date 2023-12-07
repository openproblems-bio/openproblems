#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

params_file="/tmp/datasets_openproblems_v1_params.yaml"

cat > "$params_file" << 'HERE'
param_list:
  - id: allen_brain_atlas
    obs_cell_type: label
    layer_counts: counts
    dataset_name: Mouse Brain Atlas
    dataset_url: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71585
    dataset_reference: tasic2016adult
    dataset_summary: Adult mouse primary visual cortex
    dataset_description: A murine brain atlas with adjacent cell types as assumed benchmark truth, inferred from deconvolution proportion correlations using matching 10x Visium slides (see Dimitrov et al., 2022).
    dataset_organism: mus_musculus

  - id: cengen
    obs_cell_type: cell_type
    obs_batch: experiment_code
    obs_tissue: tissue
    layer_counts: counts
    dataset_name: CeNGEN
    dataset_url: https://www.cengen.org
    dataset_reference: hammarlund2018cengen
    dataset_summary: Complete Gene Expression Map of an Entire Nervous System
    dataset_description: 100k FACS-isolated C. elegans neurons from 17 experiments sequenced on 10x Genomics.
    dataset_organism: caenorhabditis_elegans

  - id: immune_cells
    obs_cell_type: final_annotation
    obs_batch: batch
    obs_tissue: tissue
    layer_counts: counts
    dataset_name: Human immune
    dataset_url: https://theislab.github.io/scib-reproducibility/dataset_immune_cell_hum.html
    dataset_reference: luecken2022benchmarking
    dataset_summary: Human immune cells dataset from the scIB benchmarks
    dataset_description: Human immune cells from peripheral blood and bone marrow taken from 5 datasets comprising 10 batches across technologies (10X, Smart-seq2).
    dataset_organism: homo_sapiens

  - id: mouse_blood_olsson_labelled
    obs_cell_type: celltype
    layer_counts: counts
    dataset_name: Mouse myeloid
    dataset_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70245
    dataset_reference: olsson2016single
    dataset_summary: Myeloid lineage differentiation from mouse blood
    dataset_description: 660 FACS-isolated myeloid cells from 9 experiments sequenced using C1 Fluidigm and SMARTseq in 2016 by Olsson et al.
    dataset_organism: mus_musculus

  - id: mouse_hspc_nestorowa2016
    obs_cell_type: cell_type_label
    layer_counts: counts
    dataset_name: Mouse HSPC
    dataset_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81682
    dataset_reference: nestorowa2016single
    dataset_summary: Haematopoeitic stem and progenitor cells from mouse bone marrow
    dataset_description: 1656 hematopoietic stem and progenitor cells from mouse bone marrow. Sequenced by Smart-seq2. 
    dataset_organism: mus_musculus

  - id: pancreas
    obs_cell_type: celltype
    obs_batch: tech
    layer_counts: counts
    dataset_name: Human pancreas
    dataset_url: https://theislab.github.io/scib-reproducibility/dataset_pancreas.html
    dataset_reference: luecken2022benchmarking
    dataset_summary: Human pancreas cells dataset from the scIB benchmarks
    dataset_description: Human pancreatic islet scRNA-seq data from 6 datasets across technologies (CEL-seq, CEL-seq2, Smart-seq2, inDrop, Fluidigm C1, and SMARTER-seq). 
    dataset_organism: homo_sapiens

  # disabled as this is not working in openproblemsv1
  # - id: tabula_muris_senis_droplet_lung
  #   obs_cell_type: cell_type
  #   obs_batch: donor_id
  #   layer_counts: counts
  #   dataset_name: Tabula Muris Senis Lung
  #   dataset_url: https://tabula-muris-senis.ds.czbiohub.org
  #   dataset_reference: tabula2020single
  #   dataset_summary: Aging mouse lung cells from Tabula Muris Senis
  #   dataset_description: All lung cells from 10x profiles in Tabula Muris Senis, a 500k cell-atlas from 18 organs and tissues across the mouse lifespan.
  #   dataset_organism: mus_musculus

  - id: tenx_1k_pbmc
    layer_counts: counts
    dataset_name: 1k PBMCs
    dataset_url: https://www.10xgenomics.com/resources/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0
    dataset_reference: 10x2018pbmc
    dataset_summary: 1k peripheral blood mononuclear cells from a healthy donor
    dataset_description: 1k Peripheral Blood Mononuclear Cells (PBMCs) from a healthy donor. Sequenced on 10X v3 chemistry in November 2018 by 10X Genomics.
    dataset_organism: homo_sapiens

  - id: tenx_5k_pbmc
    layer_counts: counts
    dataset_name: 5k PBMCs
    dataset_url: https://www.10xgenomics.com/resources/datasets/5-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-with-cell-surface-proteins-v-3-chemistry-3-1-standard-3-1-0
    dataset_reference: 10x2019pbmc
    dataset_summary: 5k peripheral blood mononuclear cells from a healthy donor
    dataset_description: 5k Peripheral Blood Mononuclear Cells (PBMCs) from a healthy donor. Sequenced on 10X v3 chemistry in July 2019 by 10X Genomics.
    dataset_organism: homo_sapiens

  - id: tnbc_wu2021
    obs_cell_type: celltype_minor
    layer_counts: counts
    dataset_name: Triple-Negative Breast Cancer
    dataset_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118389
    dataset_reference: wu2021single
    dataset_summary: 1535 cells from six fresh triple-negative breast cancer tumors.
    dataset_description: 1535 cells from six TNBC donors by (Wu et al., 2021). This dataset includes cytokine activities, inferred using a multivariate linear model with cytokine-focused signatures, as assumed true cell-cell communication (Dimitrov et al., 2022).
    dataset_organism: homo_sapiens
    
  - id: zebrafish
    obs_cell_type: cell_type
    obs_batch: lab
    layer_counts: counts
    dataset_name: Zebrafish embryonic cells
    dataset_url: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112294
    dataset_reference: wagner2018single
    dataset_summary: Single-cell mRNA sequencing of zebrafish embryonic cells.
    dataset_description: 90k cells from zebrafish embryos throughout the first day of development, with and without a knockout of chordin, an important developmental gene. 
    dataset_organism: danio_rerio

normalization_methods: [log_cp10k, sqrt_cp10k, l1_sqrt]
output_dataset: '$id/dataset.h5ad'
output_meta: '$id/dataset_metadata.yaml'
output_state: '$id/state.yaml'
output_raw: force_null
output_normalized: force_null
output_pca: force_null
output_hvg: force_null
output_knn: force_null
publish_dir: s3://openproblems-nextflow/resources/datasets/openproblems_v1
HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision integration_build \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_openproblems_v1/main.nf \
  --workspace 53907369739130 \
  --compute-env 7IkB9ckC81O0dgNemcPJTD \
  --params-file "$params_file" \
  --config /tmp/nextflow.config