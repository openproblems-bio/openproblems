#!/bin/bash

# Run script for all test resources

echo "GEX2ADT"
viash run src/tasks/predict_modality/methods/novel/train/config.vsh.yaml -- \
  --input_train_mod1 resources_test/predict_modality/openproblems_neurips2021/bmmc_cite/normal/train_mod1.h5ad \
  --input_train_mod2 resources_test/predict_modality/openproblems_neurips2021/bmmc_cite/normal/train_mod2.h5ad \
  --output output/model.pt

# echo "ADT2GEX"
# viash run src/tasks/predict_modality/methods/novel/train/config.vsh.yaml -- \
#   --input_train_mod1 resources_test/predict_modality/openproblems_neurips2021/bmmc_cite/swap/train_mod1.h5ad \
#   --input_train_mod2 resources_test/predict_modality/openproblems_neurips2021/bmmc_cite/swap/train_mod2.h5ad \
#   --output output/model.pt

# echo "GEX2ATAC"
# viash run src/tasks/predict_modality/methods/novel/train/config.vsh.yaml -- \
#   --input_train_mod1 resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/normal/train_mod1.h5ad \
#   --input_train_mod2 resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/normal/train_mod2.h5ad \
#   --output output/model.pt

# echo "ATAC2GEX"
# viash run src/tasks/predict_modality/methods/novel/train/config.vsh.yaml -- \
#   --input_train_mod1 resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/swap/train_mod1.h5ad \
#   --input_train_mod2 resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/swap/train_mod2.h5ad \
#   --output output/model.pt


