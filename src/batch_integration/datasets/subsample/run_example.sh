SCRIPTPATH="$(
  cd "$(dirname "$0")" >/dev/null 2>&1 || exit
  pwd -P
)"

bin/viash run ${SCRIPTPATH}/config.vsh.yaml -- \
  --adata src/common/dataset_loader/download/resources/pancreas.h5ad \
  --label celltype \
  --batch tech \
  --output src/batch_integration/resources/data_loader_pancreas.h5ad \
  --debug true
