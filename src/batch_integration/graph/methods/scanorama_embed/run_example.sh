set -e
SCRIPTPATH="$(
  cd "$(dirname "$0")" >/dev/null 2>&1 || exit
  pwd -P
)"

bin/viash run ${SCRIPTPATH}/config.vsh.yaml -- \
  --input src/batch_integration/datasets/resources/datasets_pancreas.h5ad \
  --hvg true \
  --scaling true \
  --output src/batch_integration/resources/graph_pancreas_scanorama_embed.h5ad \
  --debug true
