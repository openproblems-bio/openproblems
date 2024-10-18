#!/bin/bash

## VIASH START
## VIASH END

cat > _viash.yaml << EOM
info:
  test_resources:
    - type: s3
      path: s3://openproblems-data/resources_test/common/pancreas
      dest: foo
EOM

echo ">> Run aws s3 sync"
"$meta_executable" \
  --input _viash.yaml \
  --output . \
  --quiet

echo ">> Check whether the right files were copied"
[ ! -f foo/dataset.h5ad ] && echo csv should have been copied && exit 1

echo ">> Test succeeded!"