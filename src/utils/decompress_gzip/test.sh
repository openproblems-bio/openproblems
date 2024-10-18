#!/bin/bash

set -e

## VIASH START
## VIASH END

echo "> Creating test file"
echo "Foo bar" > uncompressed.txt

echo "> Compressing file"
gzip uncompressed.txt -c > compressed.txt.gz

echo "> Decompressing file"
"$meta_executable" \
  --input "compressed.txt.gz" \
  --output "decompressed.txt"

echo "> Comparing files"
diff uncompressed.txt decompressed.txt

echo "> Test succeeded!"