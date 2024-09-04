#!/bin/bash

## VIASH START
## VIASH END

echo ">> Run aws s3 sync"
./$meta_functionality_name \
  --input s3://openproblems-data/resources_test/common/pancreas \
  --output foo \
  --quiet

echo ">> Check whether the right files were copied"
[ ! -f foo/dataset.h5ad ] && echo csv should have been copied && exit 1

echo ">> Test succeeded!"