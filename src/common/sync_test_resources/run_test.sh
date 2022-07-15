#!/bin/bash

## VIASH START
## VIASH END

echo ">> Run aws s3 sync"
./$meta_functionality_name \
  --input s3://openproblems-data/label_projection/pancreas \
  --output foo \
  --quiet

echo ">> Check whether the right files were copied"
[ ! -f foo/raw_data.h5ad ] && echo csv should have been copied && exit 1

echo ">> Test succeeded!"