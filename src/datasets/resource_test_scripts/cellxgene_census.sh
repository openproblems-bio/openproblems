#!/bin/bash

DATASET_DIR="resources_test/common/cellxgene_census"

[ ! -d "$DATASET_DIR" ] && mkdir -p "$DATASET_DIR"

# download cell ontology obo file
wget https://github.com/obophenotype/cell-ontology/releases/download/v2023-02-15/cl.obo -O "$DATASET_DIR/cl.obo"

# fetch dataset from cellxgene census
viash run src/datasets/loaders/query_cellxgene_census/config.vsh.yaml -- \
  --census_version 2023-07-25 \
  --obs_value_filter "is_primary_data == True and cell_type_ontology_term_id in ['CL:0000136', 'CL:1000311', 'CL:0002616'] and suspension_type == 'cell'" \
  --output "$DATASET_DIR/dataset.h5ad"
