#!/bin/bash
#
#make sure the following command has been executed
#viash_build -q 'label_projection|common'

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# run benchmark
export NXF_VER=23.04.2

cat > /tmp/params.yaml << HERE
param_list:
HERE

# for id in bmmc_cite_starter bmmc_cite_starter_swapped bmmc_multiome_starter bmmc_multiome_starter_swapped; do
for id in `ls resources_test/predict_modality/`; do
cat >> /tmp/params.yaml << HERE
  - id: $id
    input_train_mod1: resources_test/predict_modality/$id/train_mod1.h5ad
    input_train_mod2: resources_test/predict_modality/$id/train_mod2.h5ad
    input_test_mod1: resources_test/predict_modality/$id/test_mod1.h5ad
    input_test_mod2: resources_test/predict_modality/$id/test_mod2.h5ad
HERE
done

nextflow \
  run . \
  -main-script src/tasks/predict_modality/workflows/run/main.nf \
  -profile docker \
  -resume \
  -params-file /tmp/params.yaml \
  --output scores.tsv \
  --publish_dir output/predict_modality/