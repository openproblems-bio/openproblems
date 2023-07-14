#!/bin/bash
#
#make sure the following command has been executed
#viash ns build -q 'predict_modality|common' --parallel --setup cb

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e
# TODO: Download the starter datasets from the source repository
# TODO TODO: Generate the datasets from the source GEO dataset

generate_pm_test_resources () {
  DATASET_ID="$1"
  MOD_1_DATA="$2"
  MOD_2_DATA="$3"
  DATASET_DIR="$4"
  FLAGS="$5"

  if [ ! -f $MOD_1_DATA ]; then
      echo "Error! Could not find raw data"
      exit 1
  fi

  mkdir -p $DATASET_DIR


  # process_dataset
  viash run src/tasks/predict_modality/process_dataset/config.vsh.yaml -- \
    --input_rna $MOD_1_DATA \
    --input_other_mod $MOD_2_DATA \
    --output_train_mod1 $DATASET_DIR/train_mod1.h5ad \
    --output_train_mod2 $DATASET_DIR/train_mod2.h5ad \
    --output_test_mod1 $DATASET_DIR/test_mod1.h5ad \
    --output_test_mod2 $DATASET_DIR/test_mod2.h5ad $FLAGS

  # run one method
  viash run src/tasks/predict_modality/methods/knnr_py/config.vsh.yaml -- \
    --input_train_mod1 $DATASET_DIR/train_mod1.h5ad \
    --input_train_mod2 $DATASET_DIR/train_mod2.h5ad \
    --input_test_mod1 $DATASET_DIR/test_mod1.h5ad \
    --output $DATASET_DIR/prediction.h5ad

  # run one metric
  viash run src/tasks/predict_modality/metrics/mse/config.vsh.yaml -- \
    --input_prediction $DATASET_DIR/prediction.h5ad \
    --input_test_mod2 $DATASET_DIR/test_mod2.h5ad \
    --output $DATASET_DIR/score.h5ad

  # run benchmark on test data
  export NXF_VER=22.04.5

  nextflow run . \
    -main-script src/tasks/predict_modality/workflows/run/main.nf \
    -profile docker \
    -c src/wf_utils/labels_ci.config \
    --id "$DATASET_ID" \
    --input_train_mod1 $DATASET_DIR/train_mod1.h5ad \
    --input_train_mod2 $DATASET_DIR/train_mod2.h5ad \
    --input_test_mod1 $DATASET_DIR/test_mod1.h5ad \
    --input_test_mod2 $DATASET_DIR/test_mod2.h5ad \
    --output scores.tsv \
    --publish_dir $DATASET_DIR/
}

generate_pm_test_resources \
  bmmc_cite_starter \
  resources_test/common/bmmc_cite_starter/openproblems_bmmc_cite_starter.output_rna.h5ad \
  resources_test/common/bmmc_cite_starter/openproblems_bmmc_cite_starter.output_adt.h5ad \
  resources_test/predict_modality/bmmc_cite_starter \
  ""

generate_pm_test_resources \
  bmmc_cite_starter_swapped \
  resources_test/common/bmmc_cite_starter/openproblems_bmmc_cite_starter.output_adt.h5ad \
  resources_test/common/bmmc_cite_starter/openproblems_bmmc_cite_starter.output_rna.h5ad \
  resources_test/predict_modality/bmmc_cite_starter_swapped \
  "--swap true"

generate_pm_test_resources \
  bmmc_multiome_starter \
  resources_test/common/bmmc_multiome_starter/openproblems_bmmc_multiome_starter.output_rna.h5ad \
  resources_test/common/bmmc_multiome_starter/openproblems_bmmc_multiome_starter.output_atac.h5ad \
  resources_test/predict_modality/bmmc_multiome_starter \
  ""

generate_pm_test_resources \
  bmmc_multiome_starter_swapped \
  resources_test/common/bmmc_multiome_starter/openproblems_bmmc_multiome_starter.output_atac.h5ad \
  resources_test/common/bmmc_multiome_starter/openproblems_bmmc_multiome_starter.output_rna.h5ad \
  resources_test/predict_modality/bmmc_multiome_starter_swapped \
  "--swap true"