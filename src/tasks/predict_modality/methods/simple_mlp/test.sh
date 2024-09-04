#!/bin/bash

DATASET_DIR=resources_test/predict_modality/openproblems_neurips2021/bmmc_multiome/swap

viash run src/tasks/predict_modality/methods/simple_mlp/train/config.vsh.yaml -- \
    --input_train_mod1 $DATASET_DIR/train_mod1.h5ad \
    --input_train_mod2 $DATASET_DIR/train_mod2.h5ad \
    --output 'output/model'

viash run src/tasks/predict_modality/methods/simple_mlp/predict/config.vsh.yaml -- \
    --input_train_mod2 $DATASET_DIR/train_mod2.h5ad \
    --input_test_mod1 $DATASET_DIR/test_mod1.h5ad \
    --input_model 'output/model' \
    --output 'output/predictions'
