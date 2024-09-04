REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

nextflow run . \
  -main-script target/nextflow/predict_modality/methods/novel/main.nf \
  -profile docker \
  -c src/wf_utils/labels_ci.config \
  --input_train_mod1 resources/predict_modality/datasets/openproblems_neurips2021/bmmc_cite/normal/log_cp10k/train_mod1.h5ad \
  --input_train_mod2 resources/predict_modality/datasets/openproblems_neurips2021/bmmc_cite/normal/log_cp10k/train_mod2.h5ad \
  --input_test_mod1 resources/predict_modality/datasets/openproblems_neurips2021/bmmc_cite/normal/log_cp10k/test_mod1.h5ad \
  --publish_dir output/novel/nextflow
