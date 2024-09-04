#!/bin/bash

RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"
publish_dir="s3://openproblems-data/resources/spatially_variable_genes/results/${RUN_ID}"

# cat > /tmp/params.yaml << HERE
# input_states: s3://openproblems-data/resources/spatially_variable_genes/datasets/**/state.yaml
# rename_keys: 'input_dataset:output_dataset,input_solution:output_solution'
# output_state: "state.yaml"
# publish_dir: "$publish_dir"
# HERE

cat > /tmp/params.yaml << HERE
param_list:
  - id: svg_datasets_visium
    input_states: "s3://openproblems-data/resources/spatially_variable_genes/datasets/spatial_10x_visium/**/state.yaml"
    settings: '{"coord_type_moran_i": "generic", "coord_type_sepal": "grid", "max_neighs_sepal": 6}'

  - id: svg_datasets_xenium
    input_states: "s3://openproblems-data/resources/spatially_variable_genes/datasets/spatial_10x_xenium/**/state.yaml"
    settings: '{"coord_type_moran_i": "generic", "coord_type_sepal": "grid", "max_neighs_sepal": 6}'

  - id: svg_datasets_dbitseq
    input_states: "s3://openproblems-data/resources/spatially_variable_genes/datasets/spatial_dbit_seq/**/state.yaml"
    settings: '{"coord_type_moran_i": "generic", "coord_type_sepal": "grid", "max_neighs_sepal": 4}'

  - id: svg_datasets_merfish
    input_states: "s3://openproblems-data/resources/spatially_variable_genes/datasets/spatial_merfish/**/state.yaml"
    settings: '{"coord_type_moran_i": "generic", "coord_type_sepal": "grid", "max_neighs_sepal": 6}'

  - id: svg_datasets_seqfish
    input_states: "s3://openproblems-data/resources/spatially_variable_genes/datasets/spatial_seqfish/**/state.yaml"
    settings: '{"coord_type_moran_i": "generic", "coord_type_sepal": "grid", "max_neighs_sepal": 6}'

  - id: svg_datasets_slidetags
    input_states: "s3://openproblems-data/resources/spatially_variable_genes/datasets/spatial_slide_tags/**/state.yaml"
    settings: '{"coord_type_moran_i": "generic", "coord_type_sepal": "grid", "max_neighs_sepal": 6}'

  - id: svg_datasets_slideseqv2
    input_states: "s3://openproblems-data/resources/spatially_variable_genes/datasets/spatial_slideseq_v2/**/state.yaml"
    settings: '{"coord_type_moran_i": "generic", "coord_type_sepal": "grid", "max_neighs_sepal": 6}'

  - id: svg_datasets_starmap
    input_states: "s3://openproblems-data/resources/spatially_variable_genes/datasets/spatial_star_map/**/state.yaml"
    settings: '{"coord_type_moran_i": "generic", "coord_type_sepal": "grid", "max_neighs_sepal": 6}'

  - id: svg_datasets_stereoseq
    input_states: "s3://openproblems-data/resources/spatially_variable_genes/datasets/spatial_stereo_seq/**/state.yaml"
    settings: '{"coord_type_moran_i": "generic", "coord_type_sepal": "grid", "max_neighs_sepal": 4}'

rename_keys: 'input_dataset:output_dataset,input_solution:output_solution'
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision integration_build \
  --pull-latest \
  --main-script target/nextflow/spatially_variable_genes/workflows/run_benchmark/main.nf \
  --workspace 53907369739130 \
  --compute-env 1euVrtATIcRyy9Yc2ERbaZ \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config src/wf_utils/labels_tw.config \
