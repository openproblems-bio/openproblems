SCRIPTPATH="$(
  cd "$(dirname "$0")" >/dev/null 2>&1 || exit
  pwd -P
)"

bin/viash run ${SCRIPTPATH}/config.vsh.yaml -- \
  --adata ${SCRIPTPATH}/../resources/data_loader_pancreas.h5ad \
  --label celltype \
  --batch tech \
  --hvgs 100 \
  --output ${SCRIPTPATH}/../resources/datasets_pancreas.h5ad \
  --debug
