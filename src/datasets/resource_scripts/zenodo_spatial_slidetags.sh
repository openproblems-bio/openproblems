#!/bin/bash

cat > "/tmp/params.yaml" << 'HERE'
param_list:
  - id: zenodo_spatial_slidetags/human_cortex_slidetags
    input_data: "https://zenodo.org/records/12785822/files/slidetag_human_cortex.tar.gz?download=1"
    dataset_name: Slide-tags - Human Cortex
    dataset_url: "https://www.nature.com/articles/s41586-023-06837-4"
    dataset_summary: Slide-tags enables single-nucleus barcoding for multimodal spatial genomics.
    dataset_description: "A 100 mm2 region of the human prefrontal cortex from a neurotypical donor aged 78 years was profiled by Slide-tags."
    dataset_organism: Homo sapiens
    dataset_reference: russell2023slide
    spot_filter_min_genes: 200
    gene_filter_min_spots: 50
    remove_mitochondrial: true

  - id: zenodo_spatial_slidetags/human_skin_melanoma_slidetags
    input_data: "https://zenodo.org/records/12785822/files/slidetag_human_skin_melanoma.tar.gz?download=1"
    dataset_name: Slide-tags - Human Skin Melanoma
    dataset_url: "https://www.nature.com/articles/s41586-023-06837-4"
    dataset_summary: Slide-tags enables single-nucleus barcoding for multimodal spatial genomics.
    dataset_description: "A metastatic melanoma sample was profiled by Slide-tags."
    dataset_organism: Homo sapiens
    dataset_reference: russell2023slide
    spot_filter_min_genes: 200
    gene_filter_min_spots: 50
    remove_mitochondrial: true

  - id: zenodo_spatial_slidetags/human_tonsil_slidetags
    input_data: "https://zenodo.org/records/12785822/files/slidetag_human_tonsil.tar.gz?download=1"
    dataset_name: Slide-tags - Human Tonsil
    dataset_url: "https://www.nature.com/articles/s41586-023-06837-4"
    dataset_summary: Slide-tags enables single-nucleus barcoding for multimodal spatial genomics.
    dataset_description: "A human tonsil was profiled by Slide-tags."
    dataset_organism: Homo sapiens
    dataset_reference: russell2023slide
    spot_filter_min_genes: 200
    gene_filter_min_spots: 50
    remove_mitochondrial: true

  - id: zenodo_spatial_slidetags/mouse_embryo_slidetags
    input_data: "https://zenodo.org/records/12785822/files/slidetag_mouse_embryo.tar.gz?download=1"
    dataset_name: Slide-tags - Mouse Embryo
    dataset_url: "https://www.nature.com/articles/s41586-023-06837-4"
    dataset_summary: Slide-tags enables single-nucleus barcoding for multimodal spatial genomics.
    dataset_description: "Mouse embryo tonsil was profiled by Slide-tags."
    dataset_organism: Mus musculus
    dataset_reference: russell2023slide
    spot_filter_min_genes: 200
    gene_filter_min_spots: 50
    remove_mitochondrial: false

normalization_methods: [log_cp10k]
output_dataset: '$id/dataset.h5ad'
output_meta: '$id/dataset_metadata.yaml'
output_state: '$id/state.yaml'
output_raw: force_null
output_normalized: force_null
publish_dir: resources/datasets
HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
  withLabel: highmem {
    memory = '350GB'
  }
  withName: '.*publishStatesProc' {
    memory = '16GB'
    disk = '100GB'
  }
}
HERE

tw launch https://github.com/openproblems-bio/openproblems.git \
  --revision main_build \
  --pull-latest \
  --main-script target/nextflow/datasets/workflows/process_zenodo_spatial_slidetags/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file "/tmp/params.yaml" \
  --config /tmp/nextflow.config 
