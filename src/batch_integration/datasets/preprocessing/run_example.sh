SCRIPTPATH="$(
  cd "$(dirname "$0")" >/dev/null 2>&1 || exit
  pwd -P
)"

bin/viash run ${SCRIPTPATH}/config.vsh.yaml -- \
  --input src/batch_integration/resources/data_loader_pancreas.h5ad \
  --label celltype \
  --batch tech \
  --hvgs 100 \
  --output src/batch_integration/resources/datasets_pancreas.h5ad \
  --debug true
