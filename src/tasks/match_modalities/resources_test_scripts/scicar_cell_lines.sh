#!/bin/bash

viash run src/tasks/match_modalities/methods/harmonic_alignment/config.vsh.yaml -- \
  --input_mod1 resources_test/common/scicar_cell_lines/dataset_mod1.h5ad \
  --input_mod2 resources_test/common/scicar_cell_lines/dataset_mod2.h5ad \
  --output_mod1 resources_test/match_modalities/scicar_cell_lines/integrated_mod1.h5ad \
  --output_mod2 resources_test/match_modalities/scicar_cell_lines/integrated_mod2.h5ad