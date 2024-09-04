#!/bin/bash

cat > /tmp/params.yaml << 'HERE'
param_list:
  - id: svg_process_datasets_visium
    input_states: "s3://openproblems-data/resources/datasets/spatial_10x_visium/**/state.yaml"
    settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 50, "num_reference_genes": 200, "coord_type_proc": "grid"}'

  - id: svg_process_datasets_xenium
    input_states: "s3://openproblems-data/resources/datasets/spatial_10x_xenium/**/state.yaml"
    settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 50, "num_reference_genes": 100, "coord_type_proc": "grid"}'

  - id: svg_process_datasets_slidetags
    input_states: "s3://openproblems-data/resources/datasets/spatial_slide_tags/**/state.yaml"
    settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 50, "num_reference_genes": 50, "coord_type_proc": "grid"}'

  - id: svg_process_datasets_slideseqv2
    input_states: "s3://openproblems-data/resources/datasets/spatial_slideseq_v2/**/state.yaml"
    settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 10, "num_reference_genes": 10, "coord_type_proc": "generic"}'

  - id: svg_process_datasets_dbitseq
    input_states: "s3://openproblems-data/resources/datasets/spatial_dbit_seq/**/state.yaml"
    settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 50, "num_reference_genes": 200, "coord_type_proc": "generic"}'

  - id: svg_process_datasets_seqfish
    input_states: "s3://openproblems-data/resources/datasets/spatial_seqfish/**/state.yaml"
    settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 25, "num_reference_genes": 25, "coord_type_proc": "generic"}'

  - id: svg_process_datasets_starmap
    input_states: "s3://openproblems-data/resources/datasets/spatial_star_map/**/state.yaml"
    settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 25, "num_reference_genes": 25, "coord_type_proc": "generic"}'

  - id: svg_process_datasets_stereoseq
    input_states: "s3://openproblems-data/resources/datasets/spatial_stereo_seq/**/state.yaml"
    settings: '{"output_dataset": "$id/dataset.h5ad", "output_solution": "$id/solution.h5ad", "dataset_simulated_normalized": "$id/simulated_dataset.h5ad", "gp_k_sim": 500, "select_top_variable_genes_sim": 50, "num_reference_genes": 50, "coord_type_proc": "generic"}'

rename_keys: 'input:output_dataset'
output_state: "$id/state.yaml"
publish_dir: s3://openproblems-data/resources/spatially_variable_genes/datasets
HERE

# cat > /tmp/params.yaml << 'HERE'
# param_list:
#   - id: spatial_merfish/human_cortex_1
#     input: "s3://openproblems-data/resources/datasets/spatial_merfish/human_cortex_1/dataset.h5ad"
#     gp_k_sim: 500
#     select_top_variable_genes_sim: 25
#     num_reference_genes: 25
#     coord_type_proc: generic
  
#   - id: spatial_merfish/human_cortex_2
#     input: "s3://openproblems-data/resources/datasets/spatial_merfish/human_cortex_2/dataset.h5ad"
#     gp_k_sim: 500
#     select_top_variable_genes_sim: 50
#     num_reference_genes: 50
#     coord_type_proc: generic
  
#   - id: spatial_merfish/human_cortex_3
#     input: "s3://openproblems-data/resources/datasets/spatial_merfish/human_cortex_3/dataset.h5ad"
#     gp_k_sim: 500
#     select_top_variable_genes_sim: 50
#     num_reference_genes: 50
#     coord_type_proc: generic
  
#   - id: spatial_merfish/human_cortex_4
#     input: "s3://openproblems-data/resources/datasets/spatial_merfish/human_cortex_4/dataset.h5ad"
#     gp_k_sim: 500
#     select_top_variable_genes_sim: 50
#     num_reference_genes: 50
#     coord_type_proc: generic
  
#   - id: spatial_merfish/mouse_cortex
#     input: "s3://openproblems-data/resources/datasets/spatial_merfish/mouse_cortex/dataset.h5ad"
#     gp_k_sim: 500
#     select_top_variable_genes_sim: 25
#     num_reference_genes: 25
#     coord_type_proc: generic

# output_dataset: "$id/dataset.h5ad"
# output_solution: "$id/solution.h5ad"
# dataset_simulated_normalized: "$id/simulated_dataset.h5ad"
# output_state: "$id/state.yaml"
# publish_dir: s3://openproblems-data/resources/spatially_variable_genes/datasets
# HERE

cat > /tmp/nextflow.config << HERE
process {
  executor = 'awsbatch'
  withName:'.*publishStatesProc' {
    memory = '16GB'
    disk = '100GB'
  }
  withLabel:highmem {
    memory = '350GB'
  }
  withLabel:hightime { 
    time = 15.h 
  }
}
HERE

tw launch https://github.com/openproblems-bio/openproblems-v2.git \
  --revision integration_build \
  --pull-latest \
  --main-script target/nextflow/spatially_variable_genes/workflows/process_datasets/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file /tmp/params.yaml \
  --config /tmp/nextflow.config \
  --entry-name auto \
