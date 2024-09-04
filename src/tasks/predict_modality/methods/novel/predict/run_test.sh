#!/bin/bash

viash run src/tasks/predict_modality/methods/novel/predict/config.vsh.yaml -- \
    --input_train_mod2 'resources/predict_modality/datasets/openproblems_neurips2021/bmmc_cite/normal/log_cp10k/train_mod2.h5ad' \
    --input_test_mod1 'resources/predict_modality/datasets/openproblems_neurips2021/bmmc_cite/normal/log_cp10k/test_mod1.h5ad' \
    --input_model output/novel/model.pt \
    --input_transform output/novel/lsi_transform.pickle \
    --output 'output/novel/novel_test.h5ad'